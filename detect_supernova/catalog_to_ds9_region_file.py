import sys

from astropy.table import Table


def main():
    catalog_file = sys.argv[1]
    region_file = catalog_file[:catalog_file.rfind(".")] + ".reg"

    read_catalog_and_write_region_file(catalog_file, region_file)


def read_catalog_and_write_region_file(
    catalog_file, region_file, x_colname="x_peak", y_colname="y_peak", label_colname="peak_value"
):

    table = Table.read(catalog_file)

    head = "global point=circle\nimage\n"
    line_format_string = "point({0:0.2f},{1:0.2f}) # text={{{2:0.2f}}}\n"

    with open(region_file, "w") as outfile:
        outfile.write(head)
        for x, y, l in zip(table[x_colname], table[y_colname], table[label_colname]):
            outfile.write(line_format_string.format(x, y, l))


if __name__ == "__main__":
    main()
