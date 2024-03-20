identify CJ1.png
# convert CJ1.png -crop 1656x2059+1991+2 CJ1_out.png

# python crop.py <path/to/image> <crop_left> <crop_right> <crop_top> <crop_bottom> <name_after_crop>
python crop.py CJ1.png 2003 12 10 14 CJ1_out.png

xdg-open CJ1_out.png
