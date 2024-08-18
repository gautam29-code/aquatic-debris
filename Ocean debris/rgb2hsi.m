function histo_anlay = rgb2hsi(color_plane)
color_plane = im2double(color_plane);
r = color_plane(:, :, 1);
g = color_plane(:, :, 2);
b = color_plane(:, :, 3);


num = 0.5*((r - g) + (r - b));
den = sqrt((r - g).^2 + (r - b).*(g - b));
theta = acos(num./(den + eps));

H = theta;
H(b > g) = 2*pi - H(b > g);
H = H/(2*pi);

num = min(min(r, g), b);
den = r + g + b;
den(den == 0) = eps;
S = 1 - 3.* num./den;

H(S == 0) = 0;

I = (r + g + b)/3;


histo_anlay = cat(3, H, S, I);