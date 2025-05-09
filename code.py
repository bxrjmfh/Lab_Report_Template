def adaptive_gamma_correction(image_path, output_path, gamma_min=1.0, gamma_max=2.5, 
                             use_histogram=True, show_result=False):
    """
    使用自适应伽马矫正方法增强低光照图像
    
    参数:
        image_path: 输入图像路径
        output_path: 输出图像路径
        gamma_min: 最小伽马值
        gamma_max: 最大伽马值
        use_histogram: 是否使用直方图分析来确定伽马值
        show_result: 是否显示处理前后的对比图
    """
    # 读取图像
    img = cv2.imread(image_path)
    if img is None:
        print(f"无法读取图像: {image_path}")
        return None
    
    # 转换为灰度图用于计算亮度统计
    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    mean_brightness = np.mean(gray)
    
    if use_histogram:
        dark_threshold = 50  # 亮度阈值
        dark_ratio = np.sum(gray < dark_threshold) / gray.size
        # 计算亮度直方图
        hist = cv2.calcHist([gray], [0], None, [256], [0, 256])
        hist_norm = hist / hist.sum()  # 归一化直方图
        # 计算累积分布函数(CDF)
        cdf = hist_norm.cumsum()
        # 使用CDF的形状来确定伽马值
        # 如果暗区域像素较多，CDF在低亮度区域上升较快，应该使用较大的伽马值
        dark_weight = cdf[dark_threshold]  # 暗区域的累积权重
        # 根据暗区域权重和暗像素比例计算伽马值
        gamma = gamma_min + (gamma_max - gamma_min) * (dark_ratio * 0.7 + dark_weight * 0.3)
        gamma = min(gamma_max, max(gamma_min, gamma))  # 确保伽马值在指定范围内
    else:
        # 简单的基于平均亮度的自适应伽马值
        # 亮度越低，伽马值越高
        brightness_normalized = mean_brightness / 255.0
        gamma = gamma_max - (gamma_max - gamma_min) * brightness_normalized
    
    print(f"图像平均亮度: {mean_brightness:.2f}, 自适应伽马值: {gamma:.2f}")
    img_normalized = img.astype(np.float32) / 255.0
    img_gamma_corrected = np.power(img_normalized, 1.0/gamma)
    img_gamma_corrected = np.clip(img_gamma_corrected * 255.0, 0, 255).astype(np.uint8)
    cv2.imwrite(output_path, img_gamma_corrected)
    print(f"处理后的图像已保存到: {output_path}")
    
    # 显示处理前后的对比图
    if show_result:
        plt.figure(figsize=(12, 6))
        
        plt.subplot(1, 2, 1)
        plt.title('原始图像')
        plt.imshow(cv2.cvtColor(img, cv2.COLOR_BGR2RGB))
        plt.axis('off')
        
        plt.subplot(1, 2, 2)
        plt.title(f'伽马矫正后 (γ={gamma:.2f})')
        plt.imshow(cv2.cvtColor(img_gamma_corrected, cv2.COLOR_BGR2RGB))
        plt.axis('off')
        
        plt.tight_layout()
        plt.show()
    
    return img_gamma_corrected