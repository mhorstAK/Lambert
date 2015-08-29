import PIL
from PIL import Image

baseheight = 500
img = Image.open('LambertsTransferDiagram.png')
hpercent = (baseheight / float(img.size[1]))
wsize = int((float(img.size[0]) * float(hpercent)))
img = img.resize((wsize, baseheight), PIL.Image.ANTIALIAS)
img.show()