import squidpy as sq
import pandas as pd
from scipy.spatial import ConvexHull
import numpy as np
from shapely.geometry import Polygon, Point
import random
# import matplotlib.pyplot as plt


def make_a_buffer(shape, cell_type, epcam_buffer=3.5, other_buffer=2.5, r=False):
    x = shape["array_row"].values
    y = shape["array_col"].values
    points = np.column_stack((x, y))

    # if there are less than 3 points, just make a buffer of 1
    if len(points) < 3:
        # get the min and max x, y values
        min_x, min_y = np.min(points, axis=0)
        max_x, max_y = np.max(points, axis=0)

        # create a buffer of 1
        polygon = Polygon(
            [(min_x, min_y), (min_x, max_y), (max_x, max_y), (max_x, min_y)]
        )
    else:
        hull = ConvexHull(points, qhull_options="QJ")

        # Create the polygon from the points on the convex hull
        polygon = Polygon(points[hull.vertices])

    # Buffer the polygon
    # Increase buffer distance as required
    if "epcam" in cell_type.lower():
        if r:
            buffered_polygon = polygon.buffer(
                random.uniform(epcam_buffer - 1.25, epcam_buffer))
        else:
            buffered_polygon = polygon.buffer(epcam_buffer)
    else:
        buffered_polygon = polygon.buffer(other_buffer)

    # Calculate the min and max x, y values
    min_x, min_y, max_x, max_y = buffered_polygon.bounds

    # Create a grid of points within the bounds
    x_cords = np.arange(np.floor(min_x), np.ceil(max_x) + 1)
    y_cords = np.arange(np.floor(min_y), np.ceil(max_y) + 1)
    xy_cords = np.transpose(
        [np.tile(x_cords, len(y_cords)), np.repeat(y_cords, len(x_cords))])

    # Check which points are within the buffered polygon
    inside_points = []
    for xy in xy_cords:
        point = Point(xy)
        if buffered_polygon.contains(point):
            inside_points.append(xy)

    inside_points = np.array(inside_points)

    # Plot points inside buffered polygon
    # plt.scatter(inside_points[:,0], inside_points[:,1], color='g', s=2)

    # plt.show()

    return inside_points


def is_in_array(array_2d, element):
    # Ensure both array_2d and element are numpy arrays
    array_2d = np.array(array_2d)
    element = np.array(element)

    # Check if element is in array_2d
    return np.any(np.all(array_2d == element, axis=1))


def find_normal(spaceranger_path, annotation_path, library_name, directory="/Users/human/Downloads/"):
    random.seed(42)
    adata = sq.read.visium(spaceranger_path)
    adata.var_names_make_unique()
    adata.obs["Barcode"] = adata.obs.index
    # get x,y location based on barcodes
    barcodes = pd.read_csv(annotation_path)

    # look at epcam
    dense = np.squeeze(np.asarray(adata[:, "EPCAM"].X.todense()))
    adata.obs["EPCAM"] = dense

    barcodes = barcodes.merge(adata.obs, on="Barcode", how="right")
    condition = (barcodes['EPCAM'] != 0) & (barcodes.iloc[:, 1].isna())
    indices = barcodes[condition].index

    barcodes.loc[indices, barcodes.columns[1]] = "epcam" + indices.astype(str)
    normal_counter = 0

    def buffer_create(epcam_buffer, other_buffer, r=False):
        nonlocal barcodes
        nonlocal normal_counter
        buffers = []
        every_cell_type = barcodes.iloc[:, 1].unique()
        for cell_type in every_cell_type:
            if pd.isna(cell_type) or cell_type == "normal":
                continue
            cell_type_barcodes = barcodes[barcodes.iloc[:, 1] == cell_type]
            buffer = make_a_buffer(cell_type_barcodes, cell_type,
                                   epcam_buffer=epcam_buffer, other_buffer=other_buffer, r=r)
            if buffer.size > 0:
                buffers.append(buffer)
        combined = np.concatenate(buffers)
        for r in barcodes.itertuples():
            if not is_in_array(combined, [r.array_row, r.array_col]) and r.EPCAM == 0 and pd.isna(r[2]):
                # set the celltypespecific to normal
                barcodes.loc[r.Index, barcodes.columns[1]] = "normal"
                normal_counter += 1
        print(f"We have {normal_counter} normal cells! Yay!")
    # do this until we find 20 normal cells :)
    minimum_normal = 16
    if normal_counter < minimum_normal:
        buffer_create(3.5, 2.5)
    if normal_counter < minimum_normal:
        buffer_create(3.5, 2.5, r=True)
    if normal_counter < minimum_normal:
        buffer_create(3, 2.5)
    if normal_counter < minimum_normal:
        buffer_create(3, 2.5, r=True)
    if normal_counter < minimum_normal:
        buffer_create(2.5, 2.5, r=True)
    if normal_counter < minimum_normal:
        buffer_create(1.5, 2.5, r=True)
    if normal_counter < minimum_normal:
        raise Exception("omg my script is a failure!")

    barcodes = barcodes.dropna(subset=[barcodes.columns[1]])
    counts = barcodes.iloc[:, 1].value_counts()
    # if there is a cell type of count 1, then remove it
    for cellType in counts.index:
        if counts[cellType] == 1:
            barcodes = barcodes[barcodes.iloc[:, 1] != cellType]

    # remove excess columns
    barcodes = barcodes.drop(
        columns=["in_tissue", "EPCAM", "array_row", "array_col"])
    # barcodes.to_csv("./test.csv", index=False)

    barcodes.to_csv(
         f"{directory}{library_name}-with_normal_annotations.csv", index=False)


if __name__ == "__main__":
    find_normal('/Users/human/Downloads/HT397B1-S1H2Fs4U1Bp1/',
                "/Users/human/Downloads/HT397B1-S1H2.csv", "HT397B1-S1H2-normal")

