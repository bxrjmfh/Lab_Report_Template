def histogram_equalization(image_path, output_path):
    """
    对输入图像进行直方图均衡化处理
    
    参数:
        image_path: 输入图像路径
        output_path: 输出图像路径
    """
    # 读取图像
    img = cv2.imread(image_path)
    
    if img is None:
        print(f"无法读取图像: {image_path}")
        return
    
    # 将图像从BGR转换为YUV颜色空间
    img_yuv = cv2.cvtColor(img, cv2.COLOR_BGR2YUV)
    # 进行直方图均衡化
    img_yuv[:,:,0] = cv2.equalizeHist(img_yuv[:,:,0])
    img_output = cv2.cvtColor(img_yuv, cv2.COLOR_YUV2BGR)
    cv2.imwrite(output_path, img_output)
    print(f"处理后的图像已保存到: {output_path}")