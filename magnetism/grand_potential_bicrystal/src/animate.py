from PIL import Image
import imageio.v2 as imageio
from pathlib import Path
import re

folder_path = Path('../Data/images/')
images = list(folder_path.glob('structure*.pgm'))
images.sort(key=lambda x: int(re.findall(r'\d+', x.stem)[0]))

image_list = []
for img_path in images:
    img = Image.open(img_path)

    img_rescaled = img.resize((512, 512), Image.NEAREST)
    image_list.append(img_rescaled)

output_path = '../Data/images/animation.mov'
imageio.mimsave(output_path, image_list, fps=10)

print("Animation complete!")
