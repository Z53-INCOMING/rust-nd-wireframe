# A Beautiful N-Dimensional Polytope Wireframe Renderer With .OFF Support
I have put a lot of time into this renderer, and love using it to make renders for my website and [YouTube channel](https://www.youtube.com/@tessimal256). I have gained lots of intuition about higher dimensional shapes by using it, and I hope that you find it similarly useful. You don't have to credit me for use, but it is appreciated, you can credit me as tessimal, @tessimal256 on YouTube, linking to my [website](https://fifty-third-dimension.neocities.org), or of course, linking to this repo page. Below is a video showcasing the renderer, pre facet expansion update.

## Renderer Showcase, K_21 Polytopes

[![Video Thumbnail](https://img.youtube.com/vi/LEB5aXmeuao/0.jpg)](https://youtu.be/LEB5aXmeuao)

All of the text files mentioned in this usage tutorial can be found in ./src.

<font color=red>Warning!!!</font> With default settings, this will crash since there are no .offs in the directory. I'd add one but it's a whole thing with the .gitignore and it doesn't matter that much and it's really easy to fix on your end.

To use this, first setup the setup.txt file, and then run main.rs, and then refer to the controls section. Once you have a good angle, you can go to rotations.txt and if you really want to motions.txt and then press Enter to render an animation to an ./images folder. When that's done rendering, you'll hear the challenge win sound effect from [Arena4D](https://tessimal.itch.io/arena4d), and you can use FFmpreg to make a gif or mp4 from it.

The two FFmpreg commands I do are usually:<br>
`ffmpeg -framerate 20 -i ./images/%03d.png -c:v libx264 -crf 16 -pix_fmt yuv420p "./VIDEO NAME HERE.mp4"`<br>
`ffmpeg -framerate 20 -i ./images/%03d.png "./GIF NAME HERE.gif"`<br>
The first one is for mp4s, and the second one is for gifs. I reccomend always rendering to an mp4, since the gifs generated with this command are really bad. I'm bad with FFmpreg, sorry. You can convert the mp4s to a gif afterward using an online tool and they turn out a lot better.

Anyway, the key bits to edit are the framerate, and the number after `-crf` is the quality. Lower values are higher quality. I don't know if these will work without an images folder, so you should probably just make that folder to be sure. Also, you're supposed to run these commands in this folder, in case that wasn't already obvious.

Uh yeah I think that's about it for how to use this program. The rest will dive into the specifics of the text files and stuff.

## rotations.txt
- every line is a plane of rotation.
- first two numbers are the axes, last number is 0 for a camera space rotation and 1 for an object space rotation.

## motion.txt
- first line is the position of the object initially. if it's too short, don't worry, the axes after just won't be set, it won't crash.
- second line is the motion over the course of the whole animation, and again, it'll only apply to the first n axes.
- motion and starting position cannot be applied to the Z (2) axis.

## setup.txt
- first line is a path to the polytope you wish to load, local to this folder starting with "./".
- second line is the horizontal and vertical resolution of the animation in a single number, so all animations must be square.
- third line is the number of frames in the animation.
- fourth line is the minimum dimension
- fifth line is the facet expansion, 0.0 turns it off, 1.0 is visually identical, 0.5 - 0.9 is good.

## Axes
0 - right<br>
1 - up<br>
2 - forward<br>
3 - ana, orange<br>
4 - V+, green<br>
5+ - some high dimensional axis, they all get treated the same, they just fade out radially

"orange" and "green" here refer to 90 degree apart hues on the hue circle. the WV (3, 4) plane is rendered as a hue circle with different angles.

## Controls
LMB/MMB + move mouse - pan XZ YZ<br>
CTRL + pan - pan XW YW<br>
Z + pan - pan XV YV<br>
X + pan - pan XU YU<br>
C + pan - pan XT YT<br>
V + pan - pan XS YS<br>
B + pan - pan XR YR<br>
N + pan - pan XQ YQ<br>
<br>
scroll wheel - make image larger or smaller<br>
Control + scroll wheel - change perspective<br>
Shift + scroll wheel - change line thickness<br>
<br>
Q/A - move fade start forward and backward. everything behind this point doesn't lose any brightness, everything after it does until fade end<br>
W/S - move fade end forward and backward. everything after this point is invisible, everything behind it gets brighter and brighter until fade start<br>
E/D - expands/contracts the extra dimensional fade, aka it controls how far you see into what's perpendicular to XYZ<br>
R/F - controls the number of subdivisions of the edges<br>
<br>
0 - reloads setup.txt<br>
<br>
Enter - starts the animation, and starts writing frames to ./images. be warned- you can still do all the usual inputs while it's rendering, if you accidentally pan or something it will be visible in the animation.
