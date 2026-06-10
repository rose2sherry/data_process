import arcpy
import os


def process_landuse_tiles(lu_dir, shp_path, out_lu_path, reference_raster=None):
    # 1. 基础环境设置
    arcpy.env.overwriteOutput = True
    arcpy.env.workspace = lu_dir

    if arcpy.CheckExtension("Spatial") == "Available":
        arcpy.CheckOutExtension("Spatial")
    else:
        print("错误：无法检出 Spatial Analyst 扩展。")
        return

    print("正在递归读取切片列表（包含所有子文件夹）...")
    raster_list = []
    for root, dirs, files in os.walk(lu_dir):
        for file in files:
            # 忽略大小写，匹配所有以 .tif 或 .TIF 结尾的文件
            if file.lower().endswith(".tif"):
                raster_list.append(os.path.join(root, file))

    # 2. 核心修复：基于统一 WGS84 坐标系的无损边界筛选
    print("正在筛选与范围相交的切片 (采用统一 WGS84 坐标系比对)...")
    shp_desc = arcpy.Describe(shp_path)
    shp_extent = shp_desc.extent  # 假设输入的 SHP 是 WGS84
    wgs84_sr = arcpy.SpatialReference(4326)  # 准备 WGS84 空间参考

    selected_rasters = []
    for r in raster_list:
        ras_desc = arcpy.Describe(r)
        try:
            # 【核心逻辑】：将每一个 UTM 小切片的边界框，反向投影到 WGS84 全球坐标系
            ras_extent_wgs84 = ras_desc.extent.projectAs(wgs84_sr)

            # 统一在 WGS84 空间下进行相交判断，彻底避开 UTM 跨带形变报错
            if not shp_extent.disjoint(ras_extent_wgs84):
                selected_rasters.append(r)
        except Exception as e:
            print(f"警告：无法处理切片 {os.path.basename(r)} 的坐标转换 - {e}")

    if not selected_rasters:
        print("未找到任何切片与输入 SHP 范围相交。请检查数据范围是否重合！")
        arcpy.CheckInExtension("Spatial")
        return

    print(f"成功筛选出 {len(selected_rasters)} 个相交切片。")

    # 3. 切片预处理（转0 & 投影）
    # 使用 scratchFolder 彻底避开 GDB 数据库锁问题
    temp_workspace = arcpy.env.scratchFolder
    processed_rasters = []
    temp_con_files = []  # 记录Con转换的中间文件，方便后续清理
    temp_mosaic_name = "temp_combined_lu.tif"  # 提前定义镶嵌名称（已加.tif防止13字符报错）

    try:
        print(f"开始对 {len(selected_rasters)} 个切片进行预处理 (255转0 & 投影WGS84)...")
        for idx, ras_path in enumerate(selected_rasters):
            print(f"  正在处理 [{idx + 1}/{len(selected_rasters)}]: {os.path.basename(ras_path)}")

            # 步骤 A：将值为 255 的背景/无效值转为 0
            ras_obj = arcpy.sa.Raster(ras_path)
            ras_con = arcpy.sa.Con(ras_obj == 255, 0, ras_obj)

            # 【关键修复】：将内存栅格保存为硬盘上的实体 .tif 文件
            temp_con_name = f"tmp_con_{idx}.tif"
            temp_con_path = os.path.join(temp_workspace, temp_con_name)
            ras_con.save(temp_con_path)
            temp_con_files.append(temp_con_path)

            # 步骤 B：投影到 WGS84
            temp_proj_name = f"preprocessed_lu_{idx}.tif"  # 强制加.tif后缀
            temp_proj_path = os.path.join(temp_workspace, temp_proj_name)

            arcpy.management.ProjectRaster(
                in_raster=temp_con_path,  # 传入刚刚实体化的磁盘文件，解决 999999
                out_raster=temp_proj_path,
                out_coor_system=wgs84_sr,
                resampling_type="NEAREST"  # 土地利用等分类数据必须用 NEAREST
            )
            processed_rasters.append(temp_proj_path)

        # 4. 像素对齐环境设置 (如果提供了参考栅格)
        if reference_raster and arcpy.Exists(reference_raster):
            arcpy.env.snapRaster = reference_raster
            arcpy.env.cellSize = reference_raster
            print(f"已启用对齐环境，参考栅格：{os.path.basename(reference_raster)}")

        # 5. 镶嵌逻辑
        if len(processed_rasters) == 1:
            print("仅1个预处理后的切片，跳过镶嵌，直接准备提取...")
            target_raster = processed_rasters[0]
        else:
            print("正在进行切片镶嵌...")
            desc = arcpy.Describe(processed_rasters[0])
            pixel_type_map = {
                'U1': '1_BIT', 'U2': '2_BIT', 'U4': '4_BIT',
                'U8': '8_BIT_UNSIGNED', 'S8': '8_BIT_SIGNED',
                'U16': '16_BIT_UNSIGNED', 'S16': '16_BIT_SIGNED',
                'U32': '32_BIT_UNSIGNED', 'S32': '32_BIT_SIGNED',
                'F32': '32_BIT_FLOAT', 'F64': '64_BIT'
            }
            ptype = pixel_type_map.get(desc.pixelType, "8_BIT_UNSIGNED")

            arcpy.management.MosaicToNewRaster(
                input_rasters=processed_rasters,
                output_location=temp_workspace,
                raster_dataset_name_with_extension=temp_mosaic_name,
                coordinate_system_for_the_raster=wgs84_sr,
                pixel_type=ptype,
                cellsize=desc.meanCellWidth,
                number_of_bands=1,
                mosaic_method="FIRST",  # 分类数据重叠区域取 FIRST 或 LAST，严禁 MEAN
                mosaic_colormap_mode="FIRST"
            )
            target_raster = os.path.join(temp_workspace, temp_mosaic_name)

        # 6. 按掩膜提取
        print("正在执行按掩膜提取 (ExtractByMask)...")
        # 此时栅格大图和 SHP 均为 WGS84，提取边界完美契合
        out_extract = arcpy.sa.ExtractByMask(target_raster, shp_path)
        out_extract.save(out_lu_path)
        print(f"成功！结果已保存至: {out_lu_path}")

    except Exception as e:
        print(f"\n❌ 处理阶段失败: {e}")

    finally:
        # 7. 垃圾回收：一次性清理所有产生的临时 TIFF 文件
        print("正在清理临时预处理文件...")

        # 清理 Con转0 和 Project的中间切片
        for tmp_ras in processed_rasters + temp_con_files:
            if arcpy.Exists(tmp_ras):
                try:
                    arcpy.management.Delete(tmp_ras)
                except:
                    pass

        # 清理 临时拼接的镶嵌大图（如果有的话）
        temp_mosaic_path = os.path.join(temp_workspace, temp_mosaic_name)
        if arcpy.Exists(temp_mosaic_path):
            try:
                arcpy.management.Delete(temp_mosaic_path)
            except:
                pass

        arcpy.CheckInExtension("Spatial")


# ================= 使用示例 =================
if __name__ == "__main__":
    LU_FOLDER = r"F:\管网生成\landuse\china"  # 土地利用切片所在文件夹 (UTM投影)
    INPUT_SHP = r"G:\全国实验\356个有城市边界数据\北京市.shp"  # 城市边界 (WGS84投影)
    OUTPUT_TIF = r"F:\管网生成\lu适应\china\北京市_lu_wgs84.tif"  # 最终输出路径

    # 如果有用于对齐像素网格的参考栅格，填入其路径；没有则留空 None
    REF_RASTER = None

    process_landuse_tiles(LU_FOLDER, INPUT_SHP, OUTPUT_TIF, REF_RASTER)