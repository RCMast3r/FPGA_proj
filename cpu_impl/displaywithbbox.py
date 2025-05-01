import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import re
from sklearn.cluster import DBSCAN

def parse_line(line):
    matches = re.search("0\\ ([A-Fa-f0-9]{3})/([A-Fa-f0-9]{3})", line)
    if matches:
        return (True, [matches.group(1), matches.group(2)])
    else:
        return (False, None)

def build_image(current_image, line_to_parse, hit_coords):
    image = current_image 
    pixel_add = parse_line(line_to_parse)
    if pixel_add[0]:
        x = int(pixel_add[1][0], 16)
        y = int(pixel_add[1][1], 16)
        image[y, x] = 255
        hit_coords.append([x, y])
    return image

def load_bounding_boxes(path):
    boxes = []
    with open(path, "r") as f:
        for line in f:
            parts = line.strip().split(",")
            if len(parts) != 4:
                continue
            min_x, min_y, max_x, max_y = map(int, parts)
            boxes.append((min_x, min_y, max_x, max_y))
    return boxes

def main():
    parser = argparse.ArgumentParser(description="Display test data with bounding boxes")
    parser.add_argument("input", help="Path to the input file", default="in.txt")
    parser.add_argument("--bboxes", help="CSV file containing bounding boxes", required=True)
    args = parser.parse_args()

    start_img = np.zeros((512, 1024))
    curr_img = start_img
    hit_coords = []
    with open(args.input, "r") as infile:
        for line in infile.readlines():
            curr_img = build_image(curr_img, line, hit_coords)

    hit_coords = np.array(hit_coords)
    fig, ax = plt.subplots()
    ax.imshow(curr_img, cmap="gray")
    ax.axis('off')
    ax.set_title(args.input)

    # Load and draw bounding boxes
    boxes = load_bounding_boxes(args.bboxes)
    for (min_x, min_y, max_x, max_y) in boxes:
        width = max_x - min_x
        height = max_y - min_y
        rect = patches.Rectangle((min_x, min_y), width, height, linewidth=1.5, edgecolor='r', facecolor='none')
        ax.add_patch(rect)

    plt.show()

if __name__ == '__main__':
    main()
