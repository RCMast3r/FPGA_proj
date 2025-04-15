import argparse
import numpy as np
import matplotlib.pyplot as plt
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

def main():
    start_img = np.zeros((512, 1024))
    curr_img = start_img.copy()
    hit_coords = []

    parser = argparse.ArgumentParser(description="Script for displaying test data")
    parser.add_argument("input", help="Path to the input file", default="in.txt")
    args = parser.parse_args()

    with open(args.input, "r") as infile:
        for line in infile.readlines():
            curr_img = build_image(curr_img, line, hit_coords)

    hit_coords = np.array(hit_coords)
    if len(hit_coords) > 0:
        clustering = DBSCAN(eps=1,min_samples=1).fit(hit_coords)
        labels = clustering.labels_
        unique_labels = set(labels)
        centroids = []

        for label in unique_labels:
            if label == -1:
                continue  # Ignore noise
            cluster_points = hit_coords[labels == label]
            centroid = cluster_points.mean(axis=0)
            centroids.append(centroid)

        # Display the image and centroids
        plt.imshow(curr_img, cmap='gray')
        plt.axis('off')
        plt.title(args.input)

        centroids = np.array(centroids)
        if len(centroids) > 0:
            plt.scatter(centroids[:, 0], centroids[:, 1], c='red', marker='x')
            print("Cluster Centroids:")
            for idx, (cx, cy) in enumerate(centroids):
                print(f"Cluster {idx}: x={cx:.2f}, y={cy:.2f}")

        plt.show()
    else:
        print("No hits detected.")

if __name__ == '__main__':
    main()
