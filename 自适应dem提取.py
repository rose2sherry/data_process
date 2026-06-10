import arcpy
import os


def extract_dem_arcpy(dem_dir, shp_path, out_dem_path, reference_raster=None):
    arcpy.env.overwriteOutput = True
    arcpy.env.workspace = dem_dir

    if arcpy.CheckExtension("Spatial") == "Available":
        arcpy.CheckOutExtension("Spatial")
    else:
        print("错误：无法检出 Spatial Analyst 扩展。")
        return

    print("正在读取切片列表...")
    raster_list = [os.path.join(dem_dir, r) for r in arcpy.ListRasters("*", "TIF")]

    print("正在筛选与范围相交的切片...")
    shp_extent = arcpy.Describe(shp_path).extent
    selected_rasters = [r for r in raster_list if not shp_extent.disjoint(arcpy.Describe(r).extent)]

    if not selected_rasters:
        print("未找到任何切片与输入 SHP 范围相交。")
        return

    # 关键环境设置：如果是深度学习数据准备，建议强制设置对齐环境
    if reference_raster and arcpy.Exists(reference_raster):
        arcpy.env.snapRaster = reference_raster
        arcpy.env.cellSize = reference_raster
        print(f"已启用对齐环境，参考栅格：{os.path.basename(reference_raster)}")

    # 情况 A：如果只涉及到 1 个切片，直接进行提取，不调用拼接工具
    if len(selected_rasters) == 1:
        print("仅找到 1 个相关切片，跳过拼接步骤，直接准备提取...")
        target_raster = selected_rasters[0]

    # 情况 B：涉及到多个切片，执行拼接逻辑
    else:
        print(f"找到 {len(selected_rasters)} 个相关切片，正在进行切片拼接...")
        desc = arcpy.Describe(selected_rasters[0])

        # 修复 ERROR 000800：建立 pixelType 到工具关键字的映射
        pixel_type_map = {
            'U1': '1_BIT', 'U2': '2_BIT', 'U4': '4_BIT',
            'U8': '8_BIT_UNSIGNED', 'S8': '8_BIT_SIGNED',
            'U16': '16_BIT_UNSIGNED', 'S16': '16_BIT_SIGNED',
            'U32': '32_BIT_UNSIGNED', 'S32': '32_BIT_SIGNED',
            'F32': '32_BIT_FLOAT', 'F64': '64_BIT'
        }
        ptype = pixel_type_map.get(desc.pixelType, "32_BIT_FLOAT")  # 默认 32位浮点

        # 添加 .tif 后缀，明确告诉 ArcPy 不要使用 Grid 格式
        temp_mosaic_name = "temp_combined_dem.tif"
        # 输出到临时文件夹而不是 GDB 数据库
        temp_workspace = arcpy.env.scratchFolder

        arcpy.management.MosaicToNewRaster(
            input_rasters=selected_rasters,
            output_location=temp_workspace,
            raster_dataset_name_with_extension=temp_mosaic_name,
            coordinate_system_for_the_raster=desc.spatialReference,
            pixel_type=ptype,  # 使用转换后的正确关键字
            cellsize=desc.meanCellWidth,
            number_of_bands=1,
            mosaic_method="MEAN",
            mosaic_colormap_mode="FIRST"
        )
        target_raster = os.path.join(temp_workspace, temp_mosaic_name)

    print("正在执行按掩膜提取 (ExtractByMask)...")
    try:
        out_extract = arcpy.sa.ExtractByMask(target_raster, shp_path)
        out_extract.save(out_dem_path)
        print(f"成功！结果已保存至: {out_dem_path}")
    except Exception as e:
        print(f"裁剪阶段失败: {e}")
    finally:
        arcpy.CheckInExtension("Spatial")


# ================= 使用示例 =================
if __name__ == "__main__":
    DEM_FOLDER = r"F:\管网生成\DEM"
    INPUT_SHP = r"G:\深度学习生成管网\codex\管网数据处理\boundary_shp\boundary.shp"
    OUTPUT_TIF = r"G:\深度学习生成管网\codex\管网数据处理\data\dem\dem.tif"

    # 如果你有做好的管网标签或建筑密度图，填入此路径可实现像素对齐
    REF_RASTER = None

    extract_dem_arcpy(DEM_FOLDER, INPUT_SHP, OUTPUT_TIF, REF_RASTER)
