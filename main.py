# -*- coding: utf-8 -*-
import arcpy
import os
import numpy as np
from arcpy.sa import *


def fast_calc_building_density(shp_path, dem_ras, out_path, sub_res_factor=10):
    """建筑密度计算 (高精度栅格聚合)"""
    print(f"      正在执行高精度建筑密度计算 (细分倍数: {sub_res_factor}x)...")

    target_cell_size = float(dem_ras.meanCellWidth)
    sub_cell_size = target_cell_size / sub_res_factor

    # --- 1：子像素阶段必须捕捉到 DEM，保证网格完全嵌套
    arcpy.env.extent = dem_ras.extent
    arcpy.env.outputCoordinateSystem = dem_ras.spatialReference
    arcpy.env.cellSize = sub_cell_size
    arcpy.env.snapRaster = dem_ras  # 锁定基准，不能设为 None

    # 栅格化
    temp_bldg_sub = "memory/temp_bldg_sub"
    arcpy.conversion.PolygonToRaster(
        in_features=shp_path,
        value_field="FID",
        out_rasterdataset=temp_bldg_sub,
        cell_assignment="CELL_CENTER",
        cellsize=sub_cell_size
    )

    # 二值化
    bldg_binary_sub = Con(IsNull(Raster(temp_bldg_sub)), 0.0, 1.0)

    # 降采样聚合
    density_aggregated = Aggregate(
        in_raster=bldg_binary_sub,
        cell_factor=sub_res_factor,
        aggregation_type="MEAN"
    )

    arcpy.env.cellSize = dem_ras
    arcpy.env.snapRaster = dem_ras
    arcpy.env.extent = dem_ras.extent
    arcpy.env.mask = dem_ras

    # 用 ExtractByMask 替代 Con，防止地图代数网格溢出崩溃
    final_density = ExtractByMask(density_aggregated, dem_ras)
    final_density.save(out_path)
    print(f"      ✅ 高精度建筑密度已保存: {out_path}")


def fast_calc_road_density(shp_path, dem_ras, out_path, sub_res_factor=5):
    """极速道路密度计算 (高精度栅格聚合)"""
    print(f"      正在执行极速道路密度计算 (矩阵降维聚合，细分倍数: {sub_res_factor}x)...")

    target_cell_size = float(dem_ras.meanCellWidth)
    sub_cell_size = target_cell_size / sub_res_factor

    # 锁定 SnapRaster 保证完美对齐
    arcpy.env.extent = dem_ras.extent
    arcpy.env.outputCoordinateSystem = dem_ras.spatialReference
    arcpy.env.cellSize = sub_cell_size
    arcpy.env.snapRaster = dem_ras  # 锁定基准

    temp_road_sub = "memory/temp_road_sub"
    desc = arcpy.Describe(shp_path)
    oid_field = desc.OIDFieldName

    arcpy.conversion.PolylineToRaster(
        in_features=shp_path,
        value_field=oid_field,
        out_rasterdataset=temp_road_sub,
        cell_assignment="MAXIMUM_LENGTH",
        cellsize=sub_cell_size
    )

    road_binary_sub = Con(IsNull(Raster(temp_road_sub)), 0.0, 1.0)

    road_count_aggregated = Aggregate(
        in_raster=road_binary_sub,
        cell_factor=sub_res_factor,
        aggregation_type="SUM"
    )

    approx_length = road_count_aggregated * sub_cell_size
    density_calc = approx_length / (target_cell_size * target_cell_size)

    arcpy.env.cellSize = dem_ras
    arcpy.env.snapRaster = dem_ras
    arcpy.env.extent = dem_ras.extent
    arcpy.env.mask = dem_ras

    final_density = ExtractByMask(density_calc, dem_ras)
    final_density.save(out_path)
    print(f"      ✅ 极速道路密度已保存: {out_path}")


def calculate_twi_metrics(dem_ras, out_dir, cell_w):
    """自适应离散区域的 TWI 计算逻辑"""
    print("      正在执行 TWI 及其水文过程计算 (已启用离散区域自适应模式)...")
    global_parallel = arcpy.env.parallelProcessingFactor
    arcpy.env.parallelProcessingFactor = "0"

    valid_mask = Con(~IsNull(dem_ras), 1)
    temp_poly = os.path.join(arcpy.env.workspace, "temp_valid_regions.shp")
    arcpy.conversion.RasterToPolygon(valid_mask, temp_poly, "NO_SIMPLIFY", "Value")

    # 存储强制保存后的内存路径字符串，而不是 Raster 对象
    twi_paths = []
    spi_paths = []
    sca_paths = []

    global_extent = arcpy.env.extent
    global_mask = arcpy.env.mask

    with arcpy.da.SearchCursor(temp_poly, ["OID@", "SHAPE@"]) as cursor:
        for row in cursor:
            oid, geom = row[0], row[1]
            min_area = cell_w * cell_w * 100
            if geom.area < min_area:
                print(f"        -> 跳过微小区域 OID:{oid} (面积过小)")
                continue

            print(f"        -> 正在计算离散区域 OID:{oid} ...")
            arcpy.env.extent = geom.extent
            local_dem = ExtractByMask(dem_ras, geom)

            try:
                ill_tif = Fill(in_surface_raster=local_dem, z_limit=1)
                owdir_tif = FlowDirection(in_surface_raster=ill_tif, force_flow="NORMAL", flow_direction_type="D8")
                owacc_tif = FlowAccumulation(in_flow_direction_raster=owdir_tif, data_type="FLOAT",
                                             flow_direction_type="D8")

                cell_area = cell_w * cell_w
                SCA_tif = Con(Raster(owacc_tif) == 0, 1, Raster(owacc_tif)) * cell_area / \
                          Con(Raster(owdir_tif) == 1, cell_w,
                              Con(Raster(owdir_tif) == 4, cell_w,
                                  Con(Raster(owdir_tif) == 16, cell_w,
                                      Con(Raster(owdir_tif) == 64, cell_w,
                                          Con(Raster(owdir_tif) == 2, cell_w * SquareRoot(2),
                                              Con(Raster(owdir_tif) == 8, cell_w * SquareRoot(2),
                                                  Con(Raster(owdir_tif) == 32, cell_w * SquareRoot(2),
                                                      Con(Raster(owdir_tif) == 128, cell_w * SquareRoot(2), 0)
                                                      )
                                                  )
                                              )
                                          )
                                      )
                                  )
                              )

                slope_fil_tif = Slope(in_raster=ill_tif, output_measurement="DEGREE", z_factor=1, method="PLANAR",
                                      z_unit="METER")
                local_spi = Ln(Raster(SCA_tif) * Tan(
                    Con(Raster(slope_fil_tif) <= 0, 0.00001, Raster(slope_fil_tif) * 3.1415926 / 180)))
                local_twi = Ln(Raster(SCA_tif) / Tan(
                    Con(Raster(slope_fil_tif) <= 0, 0.00001, Raster(slope_fil_tif) * 3.1415926 / 180)))

                # 【核心修复 1】：强制执行地图代数，保存为有明确命名的内存栅格
                twi_name = f"memory/twi_{oid}"
                spi_name = f"memory/spi_{oid}"
                sca_name = f"memory/sca_{oid}"

                local_twi.save(twi_name)
                local_spi.save(spi_name)
                SCA_tif.save(sca_name)

                # 将内存路径存入列表
                twi_paths.append(twi_name)
                spi_paths.append(spi_name)
                sca_paths.append(sca_name)

            except Exception as e:
                print(f"        ⚠️ 区域 OID:{oid} 计算失败，已跳过。错误原因: {str(e)}")
                continue

    arcpy.env.extent = global_extent
    arcpy.env.mask = global_mask
    arcpy.env.parallelProcessingFactor = global_parallel

    if len(twi_paths) > 0:
        print("      正在将各区域的 TWI/SPI/SCA 结果拼接到全局栅格并强制对齐矩阵...")
        sr = dem_ras.spatialReference

        # 【核心修复 2】：不要拼接进 memory，直接落盘生成临时物理文件，避开 vsimem 崩溃漏洞
        temp_twi_mos = os.path.join(out_dir, "temp_twi_mos.tif")
        temp_spi_mos = os.path.join(out_dir, "temp_spi_mos.tif")
        temp_sca_mos = os.path.join(out_dir, "temp_sca_mos.tif")

        arcpy.management.MosaicToNewRaster(twi_paths, out_dir, "temp_twi_mos.tif", sr, "32_BIT_FLOAT", cell_w, 1,
                                           "LAST")
        arcpy.management.MosaicToNewRaster(spi_paths, out_dir, "temp_spi_mos.tif", sr, "32_BIT_FLOAT", cell_w, 1,
                                           "LAST")
        arcpy.management.MosaicToNewRaster(sca_paths, out_dir, "temp_sca_mos.tif", sr, "32_BIT_FLOAT", cell_w, 1,
                                           "LAST")

        # 强制使用基准 DEM 掩膜，死死咬住基准像元和行列数，生成最终结果
        ExtractByMask(temp_twi_mos, dem_ras).save(os.path.join(out_dir, "twi.tif"))
        ExtractByMask(temp_spi_mos, dem_ras).save(os.path.join(out_dir, "SPI_.tif"))
        ExtractByMask(temp_sca_mos, dem_ras).save(os.path.join(out_dir, "SCA_.tif"))

        # 清理落盘的临时拼接文件
        arcpy.management.Delete(temp_twi_mos)
        arcpy.management.Delete(temp_spi_mos)
        arcpy.management.Delete(temp_sca_mos)

        # 清理循环产生的内存碎片
        for p in twi_paths + spi_paths + sca_paths:
            try:
                arcpy.management.Delete(p)
            except:
                pass

        print("      ✅ 离散区域的 TWI 及衍生指标计算、拼接与严格对齐完毕！")
    else:
        print("      ❌ 未生成任何有效的 TWI 结果，请检查输入的 DEM。")


# ==========================================
# ---main
# ==========================================
def calculate_with_arcpy(dem_path, ndvi_path, pop_path, water_shp, waterway_shp, landuse_path, building_shp, road_shp,
                         buffer_dist, out_dir):
    if arcpy.CheckExtension("Spatial") == "Available":
        arcpy.CheckOutExtension("Spatial")
    else:
        print("❌ 未获得 Spatial Analyst 许可")
        return

    arcpy.env.overwriteOutput = True
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    arcpy.env.workspace = out_dir

    try:
        print("🚀 开始处理流程...")

        # --- 0. 坐标系检查与自动 UTM 投影
        desc = arcpy.Describe(dem_path)
        sr = desc.spatialReference
        processed_dem_path = dem_path

        if sr.type == "Geographic":
            print("⚠️ 检测到 DEM 为地理坐标系，正在自动投影至 UTM 以保证距离与密度计算的准确性...")
            extent = desc.extent
            lon = (extent.XMin + extent.XMax) / 2
            lat = (extent.YMin + extent.YMax) / 2
            utm_zone = int((lon + 180) / 6) + 1
            is_south = "S" if lat < 0 else "N"
            target_sr = arcpy.SpatialReference(f"WGS 1984 UTM Zone {utm_zone}{is_south}")

            processed_dem_path = os.path.join(out_dir, "elevation.tif")
            arcpy.management.ProjectRaster(dem_path, processed_dem_path, target_sr)
            print(f"   ✅ DEM 已重投影至: WGS 1984 UTM Zone {utm_zone}{is_south}")
        else:
            print(f"   ✅ DEM 已是投影坐标系 ({sr.name})，无需重投影。")

        dem_ras = Raster(processed_dem_path)
        base_sr = dem_ras.spatialReference
        cell_w = dem_ras.meanCellWidth

        arcpy.env.outputCoordinateSystem = base_sr
        arcpy.env.snapRaster = dem_ras
        arcpy.env.extent = dem_ras.extent
        arcpy.env.cellSize = dem_ras
        arcpy.env.mask = dem_ras

        print("   ✅ 全局环境基准（范围、像元、捕捉栅格、坐标系）已完全锁定至 DEM。")

        # --- 1. 地形基础指标计算
        print("正在计算基础地形指标 (Slope/Aspect/Curvature/TRI)...")
        Con(~IsNull(dem_ras), Slope(dem_ras, "DEGREE")).save(os.path.join(out_dir, "slope.tif"))
        Con(~IsNull(dem_ras), Aspect(dem_ras)).save(os.path.join(out_dir, "aspect.tif"))
        Con(~IsNull(dem_ras), Curvature(dem_ras)).save(os.path.join(out_dir, "curvature.tif"))
        Con(~IsNull(dem_ras), FocalStatistics(dem_ras, NbrRectangle(3, 3, "CELL"), "RANGE")).save(
            os.path.join(out_dir, "tri.tif"))

        # --- 2. TWI 相关水文指标计算
        print("正在执行 TWI 等水文指标计算模块...")
        calculate_twi_metrics(dem_ras, out_dir, cell_w)

        # --- 3. FVC 对齐 (直接投影重采样 + 异常值过滤)
        # print("正在对齐 FVC 数据 (采用 BILINEAR 重采样以保证连续性)...")
        # arcpy.management.ProjectRaster(fvc_path, "memory/fvc_proj", base_sr, "BILINEAR", cell_w)
        # fvc_proj_ras = Raster("memory/fvc_proj")
        #
        # # 【过滤插值产生的极值污染】
        # # 如果你的 FVC 数据范围是 0 到 100，请将下面的 1 更改为 100
        # valid_fvc = Con((fvc_proj_ras >= 0) & (fvc_proj_ras <= 1), fvc_proj_ras)
        #
        # # 用 DEM 有效区域进行最后掩膜并保存
        # Con(~IsNull(dem_ras), valid_fvc).save(os.path.join(out_dir, "fvc.tif"))
        # print("   ✅ FVC 对齐完毕，已剔除异常大值。")

        # --- 3. FVC 计算
        print("正在处理 NDVI 并计算 FVC...")
        arcpy.management.ProjectRaster(ndvi_path, "memory/ndvi_proj", base_sr, "BILINEAR", cell_w)
        ndvi_ras = Con(~IsNull(dem_ras), Raster("memory/ndvi_proj"))

        ndvi_arr = arcpy.RasterToNumPyArray(ndvi_ras, nodata_to_value=-9999)
        valid_ndvi = ndvi_arr[(ndvi_arr > -1) & (ndvi_arr <= 1)]
        if valid_ndvi.size > 0:
            v_min = np.percentile(valid_ndvi, 5)
            v_max = np.percentile(valid_ndvi, 95)
            fvc = (ndvi_ras - v_min) / (v_max - v_min)
            out_fvc = Con(~IsNull(dem_ras), Con(fvc < 0, 0, Con(fvc > 1, 1, fvc)))
            out_fvc.save(os.path.join(out_dir, "fvc.tif"))
        else:
            print("⚠️ NDVI 数据在研究区内无效，跳过 FVC 提取。")

        # --- 4. 人口与土地利用对齐
        print("正在对齐人口数据 (采用 BILINEAR 重采样)...")
        arcpy.management.ProjectRaster(pop_path, "memory/pop_proj", base_sr, "BILINEAR", cell_w)
        Con(~IsNull(dem_ras), Raster("memory/pop_proj")).save(os.path.join(out_dir, "pop_den.tif"))

        print("正在对齐土地利用数据 (采用 NEAREST 最邻近重采样以保持分类值不被破坏)...")
        arcpy.management.ProjectRaster(landuse_path, "memory/lu_proj", base_sr, "NEAREST", cell_w)
        Con(~IsNull(dem_ras), Int(Raster("memory/lu_proj"))).save(os.path.join(out_dir, "lulc_f.tif"))

        # --- 5. 距水体距离计算 (面状水体 + 线状水系)
        print("正在计算距水体距离 (合并计算面状与线状水体)...")
        try:
            buf_val = float(buffer_dist.split(" ")[0])
        except:
            buf_val = 5000.0

        dem_ext = dem_ras.extent
        expanded_extent = arcpy.Extent(
            dem_ext.XMin - buf_val, dem_ext.YMin - buf_val,
            dem_ext.XMax + buf_val, dem_ext.YMax + buf_val,
            spatial_reference=base_sr
        )
        bounds_shp = os.path.join(out_dir, "temp_bounds_utm.shp")
        if arcpy.Exists(bounds_shp):
            arcpy.management.Delete(bounds_shp)
        arcpy.management.CreateFeatureclass(out_dir, "temp_bounds_utm.shp", "POLYGON", spatial_reference=base_sr)
        with arcpy.da.InsertCursor(bounds_shp, ["SHAPE@"]) as cursor:
            cursor.insertRow([expanded_extent.polygon])

        # 分别创建面状和线状图层
        if arcpy.Exists("water_poly_lyr"): arcpy.management.Delete("water_poly_lyr")
        if arcpy.Exists("water_line_lyr"): arcpy.management.Delete("water_line_lyr")
        arcpy.management.MakeFeatureLayer(water_shp, "water_poly_lyr")
        arcpy.management.MakeFeatureLayer(waterway_shp, "water_line_lyr")

        print("      正在执行面状与线状水体的动态空间相交查询...")
        arcpy.management.SelectLayerByLocation("water_poly_lyr", "INTERSECT", bounds_shp)
        arcpy.management.SelectLayerByLocation("water_line_lyr", "INTERSECT", bounds_shp)

        global_mask = arcpy.env.mask
        arcpy.env.mask = None
        arcpy.env.extent = expanded_extent
        global_parallel = arcpy.env.parallelProcessingFactor
        arcpy.env.parallelProcessingFactor = "0"

        try:
            poly_count = int(arcpy.management.GetCount("water_poly_lyr")[0])
            line_count = int(arcpy.management.GetCount("water_line_lyr")[0])

            dist_rasters = []

            # 1. 计算面状水体距离
            if poly_count > 0:
                print(f"      成功选中 {poly_count} 个面状水体要素...")
                dist_rasters.append(EucDistance("water_poly_lyr", cell_size=cell_w))
            else:
                print("      ⚠️ 警告: 扩展范围内未找到面状水体！")

            # 2. 计算线状水体距离
            if line_count > 0:
                print(f"      成功选中 {line_count} 个线状水体(waterways)要素...")
                dist_rasters.append(EucDistance("water_line_lyr", cell_size=cell_w))
            else:
                print("      ⚠️ 警告: 扩展范围内未找到线状水体！")

            arcpy.env.extent = dem_ras.extent
            arcpy.env.mask = global_mask

            # 3. 合并逻辑：取最小值
            if len(dist_rasters) == 2:
                print("      正在合并面状与线状水体距离 (取最小值)...")
                # 使用 CellStatistics 取两者的最小值
                combined_dist = CellStatistics(dist_rasters, "MINIMUM", "DATA")
                final_dist = Con(~IsNull(dem_ras), combined_dist)
            elif len(dist_rasters) == 1:
                final_dist = Con(~IsNull(dem_ras), dist_rasters[0])
            else:
                print("      ❌ 范围内完全没有任何水体要素，输出空栅格。")
                final_dist = Con(~IsNull(dem_ras), 99999)  # 兜底值

            final_dist.save(os.path.join(out_dir, "distowater.tif"))
            print("      ✅ 距水体(面+线)最终距离已计算并成功保存！")

        except Exception as e:
            print(f"      ❌ 距水体距离计算失败: {str(e)}")
        finally:
            if arcpy.Exists("water_poly_lyr"): arcpy.management.Delete("water_poly_lyr")
            if arcpy.Exists("water_line_lyr"): arcpy.management.Delete("water_line_lyr")
            if arcpy.Exists(bounds_shp):
                try:
                    arcpy.management.Delete(bounds_shp)
                except:
                    pass
            arcpy.env.extent = dem_ras.extent
            arcpy.env.mask = global_mask
            arcpy.env.parallelProcessingFactor = global_parallel

        # --- 6. 密度分析
        print("正在执行密度分析...")
        fast_calc_building_density(building_shp, dem_ras, os.path.join(out_dir, "build_den.tif"))
        fast_calc_road_density(road_shp, dem_ras, os.path.join(out_dir, "road_den.tif"))

        # --- 7. 管网标签生成
        # print("正在生成管网掩膜标签...")
        # temp_p_ras = "memory/pipe_raster"
        # arcpy.conversion.PolylineToRaster(pipe_path, "FID", temp_p_ras, "MAXIMUM_LENGTH", "NONE", cell_w)
        # pipe_binary = Con(IsNull(Raster(temp_p_ras)), 0, 1)
        #
        # arcpy.env.compression = "LZW"
        # final_pipe = Con(~IsNull(dem_ras), pipe_binary)
        # final_pipe.save(os.path.join(out_dir, "pipenet.tif"))
        # arcpy.env.compression = "NONE"

        print("\n✨ 所有指标与标签任务处理成功！特征数据矩阵大小已完美锁定。")

    except Exception as e:
        print(f"\n❌ 运行时发生错误: {str(e)}")
        import traceback
        traceback.print_exc()
    finally:
        # 清理内存工作空间，防止内存泄漏
        arcpy.management.Delete("memory")
        arcpy.CheckInExtension("Spatial")


if __name__ == "__main__":
    config = {
        "dem": r"F:\管网生成\dem适应\china\北京市_dem.tif",# 需更换
        # "fvc": r"F:\管网生成\fvc\xyt_max_fvc_2020.tif",
        "ndvi": r"F:\管网生成\NDVI\NDVImax2020_float_final.tif",
        "pop": r"F:\管网生成\pop\chn_ppp_2020_UNadj.tif",
        "water": r"F:\管网生成\水体矢量\中国\gis_osm_water_a_free_1.shp",
        "waterway": r"F:\管网生成\水体矢量\中国\gis_osm_waterways_free_1.shp",
        "landuse": r"F:\管网生成\lu适应\china\北京市_lu_wgs84.tif", # 需更换
        # "pipe": r"F:\管网生成\管网\西雅图\SPU_DWW_Mainlines_Pe_wgs.shp",
        "building": r"F:\管网生成\gee_bulid\china\beijing.shp", # 需更换
        "road": r"F:\管网生成\水体矢量\中国\gis_osm_roads_free_1.shp",
        "buffer": "5000 Meters",
        "out": r"F:\管网生成\样本生成\china\beijing" # 需更换输出路径名称
    }

    calculate_with_arcpy(
        config["dem"], config["ndvi"], config["pop"],
        config["water"], config["waterway"], config["landuse"],
        config["building"], config["road"],
        config["buffer"], config["out"]
    )