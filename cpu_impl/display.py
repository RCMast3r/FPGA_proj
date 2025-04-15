import argparse
import numpy as np
import matplotlib.pyplot as plt
import re
# this python script displays the data in grid format with the hits / misses 

def parse_line(line):

    matches = re.search("0\\ ([A-Fa-f0-9]{3})/([A-Fa-f0-9]{3})", line)

    if matches:
        return (True, [matches.group(1), matches.group(2)])
    else:
        return (False, None)

def build_image(current_image, line_to_parse):

    image = current_image 
    pixel_add = parse_line(line_to_parse)
    if(pixel_add[0]):
        image[int(pixel_add[1][1], 16), int(pixel_add[1][0], 16)] = 255

    return image

def main():

    start_img = np.zeros((512,1024))
    curr_img = start_img
    parser = argparse.ArgumentParser(description="Script for displaying test data")
    parser.add_argument("input",  help="Path to the input file", default="in.txt")

    args = parser.parse_args()

    with open(args.input, "r") as infile:
        for line in infile.readlines():
            curr_img = build_image(curr_img, line)

    plt.imshow(curr_img)
    plt.axis('off')  # Hide axis ticks
    plt.title(args.input)
    plt.show()


if __name__ =='__main__':
    main()
