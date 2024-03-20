identify CJ2.png
# python crop.py <path/to/image> <crop_left> <crop_right> <crop_top> <crop_bottom> <name_after_crop>
python crop.py CJ2.png 6 883 510 460 CJ2_e.png
python crop.py CJ2.png 879 10 510 460 CJ2_f.png
python crop.py CJ2.png 6 883 958 12 CJ2_g.png

convert -append CJ2_e.png CJ2_f.png CJ2_g.png CJ2_out.png

xdg-open CJ2_out.png
