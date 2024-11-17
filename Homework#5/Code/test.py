import numpy as np
import cv2

def bilinear_interpolation(img, x, y):
    h, w = img.shape
    x_int = int(np.floor(x))
    y_int = int(np.floor(y))

    x_int = min(max(x_int, 0), w - 2)
    y_int = min(max(y_int, 0), h - 2)
    
    h1 = y - y_int
    h2 = 1 - h1
    w1 = x - x_int
    w2 = 1 - w1
    
    A = img[y_int, x_int]
    B = img[y_int, x_int + 1]
    C = img[y_int + 1, x_int]
    D = img[y_int + 1, x_int + 1]
    
    P = (w2 * (h2 * A + h1 * C) + w1 * (h2 * B + h1 * D))
    
    return np.clip(P, 0, 255)

def resample_image(image_path, target_width, target_height):
    image = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)

    original_height, original_width = image.shape
    print(f"Original resolution: {original_width}x{original_height}")
    
    scale_x = original_width / target_width
    scale_y = original_height / target_height
    
    resampled_image = np.zeros((target_height, target_width), dtype=np.uint8)
    
    for i in range(target_height):
        for j in range(target_width):
            x = j * scale_x
            y = i * scale_y
            resampled_image[i, j] = bilinear_interpolation(image, x, y)
    
    resized_original = cv2.resize(image, (target_width, target_height), interpolation=cv2.INTER_LINEAR)
    combined_image = cv2.hconcat([resized_original, resampled_image])
    
    cv2.imshow("원본 이미지와 리샘플링된 이미지 (이차원 보간)", combined_image)
    cv2.waitKey(0)
    cv2.destroyAllWindows()

def main():
    target_width = int(input("Enter target width (M'): "))
    target_height = int(input("Enter target height (N'): "))

    image_path = "bin.jpg"
    resample_image(image_path, target_width, target_height)

if __name__ == "__main__":
    main()
