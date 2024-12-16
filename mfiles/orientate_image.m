function Image = orientate_image(Image,orientation)

switch orientation
    case 1
    case 2
        Image = rot90(Image,1);
    case 3
        Image = rot90(Image,2);
    case 4
        Image = rot90(Image,3);
    case 5
        Image = flipud(Image);
    case 6
        Image = flipud(rot90(Image,1));
    case 7
        Image = flipud(rot90(Image,2));
    case 8
        Image = flipud(rot90(Image,3));
end