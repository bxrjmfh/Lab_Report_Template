def clahe_yuv_enhancement(image_path, output_path, clip_limit=2.0, tile_grid_size=(4, 4)):
    """
    在YUV色彩空间中使用CLAHE方法进行图像增强
    
    参数:
        image_path: 输入图像路径
        output_path: 输出图像路径
        clip_limit: CLAHE的对比度限制参数，默认为2.0
        tile_grid_size: CLAHE的网格大小，默认为8x8
    """
    # 读取图像
    img = cv2.imread(image_path)
    
    if img is None:
        print(f"无法读取图像: {image_path}")
        return
    
    # 将图像从BGR转换为YUV颜色空间
    img_yuv = cv2.cvtColor(img, cv2.COLOR_BGR2YUV)
    y, u, v = cv2.split(img_yuv)
    
    # 创建CLAHE对象
    clahe = cv2.createCLAHE(clipLimit=clip_limit, tileGridSize=tile_grid_size)
    
    # 对Y通道应用CLAHE
    y_clahe = clahe.apply(y) 
    # 合并通道
    img_yuv_enhanced = cv2.merge([y_clahe, u, v])
    img_enhanced = cv2.cvtColor(img_yuv_enhanced, cv2.COLOR_YUV2BGR)
    
    # 保存结果图像
    cv2.imwrite(output_path, img_enhanced)
    print(f"处理后的图像已保存到: {output_path}")
    
    return img_enhanced