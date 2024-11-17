import numpy as np
import cv2

def bilinear_interpolation(img, x, y):
    h, w, c = img.shape
    x_int = int(np.floor(x))
    y_int = int(np.floor(y))
    
    x_int = min(max(x_int, 0), w - 2)
    y_int = min(max(y_int, 0), h - 2)
    
    h1 = y - y_int
    h2 = 1 - h1
    w1 = x - x_int
    w2 = 1 - w1
    
    #deal with colour values
    result = np.zeros(3)
    for channel in range(c):
        A = img[y_int, x_int, channel]
        B = img[y_int, x_int + 1, channel]
        C = img[y_int + 1, x_int, channel]
        D = img[y_int + 1, x_int + 1, channel]

        result[channel] = (w2 * (h2 * A + h1 * C) + w1 * (h2 * B + h1 * D))
    
    return np.clip(result, 0, 255).astype(np.uint8)

def resample_image(image_path, target_width, target_height):
    image = cv2.imread(image_path, cv2.IMREAD_COLOR)

    original_height, original_width, _ = image.shape
    print(f"original resolution: {original_width}*{original_height}")
    
    scale_x = original_width / target_width
    scale_y = original_height / target_height

    resampled_image = np.zeros((target_height, target_width, 3), dtype=np.uint8)
    
    for i in range(target_height):
        for j in range(target_width):
            x = j * scale_x
            y = i * scale_y
            resampled_image[i, j] = bilinear_interpolation(image, x, y)
    
    cv2.imshow("Original Image", image)
    cv2.imshow("Resampled Image", resampled_image)
    
    cv2.waitKey(0)
    cv2.destroyAllWindows()

def main():
    target_width = int(input("Enter target width (M'): "))
    target_height = int(input("Enter target height (N'): "))
    
    image_path = "pik.jpg"
    resample_image(image_path, target_width, target_height)

if __name__ == "__main__":
    main()
