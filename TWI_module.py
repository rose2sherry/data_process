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
    """
    自适应离散区域的 TWI 计算逻辑 (支持多块不连通的 DEM)
    """
    print("      正在执行 TWI 及其水文过程计算 (已启用离散区域自适应模式)...")

    # --- 1. 识别并提取离散的有效区域 ---
    valid_mask = Con(~IsNull(dem_ras), 1)
    temp_poly = os.path.join(arcpy.env.workspace, "temp_valid_regions.shp")

    # 将栅格有效区域转为多边形（不简化边界，保持栅格边缘）
    arcpy.conversion.RasterToPolygon(valid_mask, temp_poly, "NO_SIMPLIFY", "Value")

    twi_rasters = []
    spi_rasters = []
    sca_rasters = []

    # 备份全局环境
    global_extent = arcpy.env.extent
    global_mask = arcpy.env.mask

    # --- 2. 遍历每一个离散块进行计算
    with arcpy.da.SearchCursor(temp_poly, ["OID@", "SHAPE@"]) as cursor:
        for row in cursor:
            oid = row[0]
            geom = row[1]

            # 过滤掉太小的离散“噪点”区域（例如小于 100 个像元的碎块）
            # 这种微小斑块做水文分析毫无物理意义，且极易导致 Fill 工具报错
            min_area = cell_w * cell_w * 100
            if geom.area < min_area:
                print(f"        -> 跳过微小区域 OID:{oid} (面积过小)")
                continue

            print(f"        -> 正在计算离散区域 OID:{oid} ...")

            # 锁定当前碎块的局部环境
            arcpy.env.extent = geom.extent

            # 提取局部 DEM
            local_dem = ExtractByMask(dem_ras, geom)

            try:
                # -- 局部水文计算开始 --
                ill_tif = Fill(in_surface_raster=local_dem, z_limit=1)
                owdir_tif = FlowDirection(in_surface_raster=ill_tif, force_flow="NORMAL", flow_direction_type="D8")
                owacc_tif = FlowAccumulation(in_flow_direction_raster=owdir_tif, data_type="FLOAT",
                                             flow_direction_type="D8")

                # SCA
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

                # 局部 SPI 和 TWI
                local_spi = Ln(Raster(SCA_tif) * Tan(
                    Con(Raster(slope_fil_tif) <= 0, 0.00001, Raster(slope_fil_tif) * 3.1415926 / 180)))
                local_twi = Ln(Raster(SCA_tif) / Tan(
                    Con(Raster(slope_fil_tif) <= 0, 0.00001, Raster(slope_fil_tif) * 3.1415926 / 180)))

                # 存入内存，供后续拼接
                twi_rasters.append(local_twi)
                spi_rasters.append(local_spi)
                sca_rasters.append(SCA_tif)
                # -- 局部水文计算结束 --

            except Exception as e:
                print(f"        ⚠️ 区域 OID:{oid} 计算失败，已跳过。错误原因: {str(e)}")
                continue

    arcpy.env.extent = global_extent
    arcpy.env.mask = global_mask

    # --- 3. 镶嵌 (Mosaic) 所有局部的计算结果
    if len(twi_rasters) > 0:
        print("      正在将各区域的 TWI/SPI/SCA 结果拼接到全局栅格...")

        # 定义拼接的公共参数
        sr = dem_ras.spatialReference

        # 拼接 TWI
        arcpy.management.MosaicToNewRaster(
            input_rasters=twi_rasters,
            output_location=out_dir,
            raster_dataset_name_with_extension="twi.tif",
            coordinate_system_for_the_raster=sr,
            pixel_type="32_BIT_FLOAT",
            cellsize=cell_w,
            number_of_bands=1,
            mosaic_method="LAST"
        )

        # 拼接 SPI
        arcpy.management.MosaicToNewRaster(
            input_rasters=spi_rasters,
            output_location=out_dir,
            raster_dataset_name_with_extension="SPI_.tif",
            coordinate_system_for_the_raster=sr,
            pixel_type="32_BIT_FLOAT",
            cellsize=cell_w,
            number_of_bands=1,
            mosaic_method="LAST"
        )

        # 拼接 SCA
        arcpy.management.MosaicToNewRaster(
            input_rasters=sca_rasters,
            output_location=out_dir,
            raster_dataset_name_with_extension="SCA_.tif",
            coordinate_system_for_the_raster=sr,
            pixel_type="32_BIT_FLOAT",
            cellsize=cell_w,
            number_of_bands=1,
            mosaic_method="LAST"
        )
        print("      ✅ 离散区域的 TWI 及衍生指标计算、拼接完毕！")
    else:
        print("      ❌ 未生成任何有效的 TWI 结果，请检查输入的 DEM。")


# ==========================================
# ---main
# ==========================================
def calculate_with_arcpy(dem_path, ndvi_path, pop_path, water_shp, landuse_path, pipe_path, building_shp, road_shp,
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

        # 将基准 DEM 实例化为 Raster
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

        # --- 4. 人口与土地利用
        print("正在重采样人口数据 (BILINEAR)...")
        arcpy.management.ProjectRaster(pop_path, "memory/pop_proj", base_sr, "BILINEAR", cell_w)
        Con(~IsNull(dem_ras), Raster("memory/pop_proj")).save(os.path.join(out_dir, "pop_den.tif"))

        print("正在重采样土地利用数据 (NEAREST 保持分类值)...")
        arcpy.management.ProjectRaster(landuse_path, "memory/lu_proj", base_sr, "NEAREST", cell_w)
        Con(~IsNull(dem_ras), Int(Raster("memory/lu_proj"))).save(os.path.join(out_dir, "lulc_f.tif"))

        # --- 5. 距水体距离计算
        print("正在计算距水体距离 (采用内存图层空间查询，规避底层验证)...")

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

        if arcpy.Exists("water_lyr"):
            arcpy.management.Delete("water_lyr")
        arcpy.management.MakeFeatureLayer(water_shp, "water_lyr")

        print("      正在执行动态空间相交查询...")
        arcpy.management.SelectLayerByLocation("water_lyr", "INTERSECT", bounds_shp)

        global_mask = arcpy.env.mask
        arcpy.env.mask = None
        arcpy.env.extent = expanded_extent

        global_parallel = arcpy.env.parallelProcessingFactor
        arcpy.env.parallelProcessingFactor = "0"

        try:
            water_count = int(arcpy.management.GetCount("water_lyr")[0])
            if water_count == 0:
                print("      ⚠️ 警告: 扩展矩形范围内未找到任何水体要素！")
            else:
                print(f"      成功选中 {water_count} 个水体要素，启动底层运算...")
                # EucDistance 会自动且仅对“被选中”的要素进行计算
                out_euc_dist = EucDistance("water_lyr", cell_size=cell_w)

                arcpy.env.extent = dem_ras.extent
                arcpy.env.mask = global_mask

                # 利用离散 DEM 有效值进行最终裁剪
                final_dist = Con(~IsNull(dem_ras), out_euc_dist)
                final_dist.save(os.path.join(out_dir, "distowater.tif"))
                print("      ✅ 距水体距离已计算并成功保存！")

        except Exception as e:
            print(f"      ❌ 距水体距离计算失败: {str(e)}")

        finally:
            # 清理虚拟图层和临时框
            if arcpy.Exists("water_lyr"):
                arcpy.management.Delete("water_lyr")
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
        fast_calc_building_density(building_shp, dem_ras, os.path.join(out_dir, "building_den.tif"))
        fast_calc_road_density(road_shp, dem_ras, os.path.join(out_dir, "road_den.tif"))

        # --- 7. 管网标签生成
        print("正在生成管网掩膜标签...")
        temp_p_ras = "memory/pipe_raster"
        # 直接利用 MAXIMUM_LENGTH 转换
        arcpy.conversion.PolylineToRaster(pipe_path, "FID", temp_p_ras, "MAXIMUM_LENGTH", "NONE", cell_w)

        # 将无管网区域补 0，并用 DEM 边界截断，生成 0/1 标签
        pipe_binary = Con(IsNull(Raster(temp_p_ras)), 0, 1)

        arcpy.env.compression = "LZW"  # 标签数据可使用 LZW 压缩降低体积
        final_pipe = Con(~IsNull(dem_ras), pipe_binary)
        final_pipe.save(os.path.join(out_dir, "pipenet.tif"))
        arcpy.env.compression = "NONE"  # 恢复默认

        print("\n✨ 所有指标与标签任务处理成功！特征数据矩阵大小已完美锁定。")

    except Exception as e:
        print(f"\n❌ 运行时发生错误: {str(e)}")
        import traceback
        traceback.print_exc()
    finally:
        arcpy.management.Delete("memory")
        arcpy.CheckInExtension("Spatial")


if __name__ == "__main__":
    config = {
        "dem": r"F:\管网生成\dem适应\hangzhou_dem.tif",
        "ndvi": r"F:\管网生成\NDVI\NDVImax2020_float_final.tif",
        "pop": r"E:\personal\dl_data\cn\chn_ppp_2020_UNadj.tif",
        "water": r"G:\1111\data_o\osm_water_all_2020.shp",
        "landuse": r"E:\personal\dl_data\cn\lucc_landsue_2020.tif",
        "pipe": r"G:\管网\sc\Link.shp",
        "building": r"F:\管网生成\gee_bulid\hangzhou.shp",
        "road": r"E:\global\全国实验\道路文件处理\osm\gis_osm_roads_all_china.shp",
        "buffer": "5000 Meters",
        "out": r"G:\1111\results\5"  #结果文件夹
    }

    calculate_with_arcpy(
        config["dem"], config["ndvi"], config["pop"],
        config["water"], config["landuse"], config["pipe"],
        config["building"], config["road"],
        config["buffer"], config["out"]
    )