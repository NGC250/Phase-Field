import matplotlib.pyplot as plt
import pandas as pd
import imageio.v2 as iio
from PIL import Image
from pathlib import Path
from joblib import Parallel, delayed
import time
from natsort import natsorted

start = time.time()

H = 256
W = 256
dt = 0.1
iteration = 100000
save_after = 500

read_dir = "../Data/status/"
write_dir = "../Data/images/"
wd_mx = write_dir + "mx/"
wd_my = write_dir + "my/"
wd_mz = write_dir + "mz/"

def image_save(i, rd):
    
    try:
        data = pd.read_csv(f"{rd}domains{i}.dat", delimiter=",")
        
        img1 = (data.iloc[:, 0]).to_numpy().reshape((H, W))
        img2 = (data.iloc[:, 1]).to_numpy().reshape((H, W))
        img3 = (data.iloc[:, 2]).to_numpy().reshape((H, W))

        plt.imshow(img1, cmap='jet', vmin=np.min(img1), vmax=np.max(img1))
        plt.colorbar()
        plt.title(f"mx domain status: t = {(i * dt):.2f}\n(dt = {dt})", pad=10)
        plt.savefig(f"{wd_mx}domains_mx{i}.png")
        plt.clf()
        plt.close(plt.gcf())

        plt.imshow(img2, cmap='jet', vmin=np.min(img2), vmax=np.max(img2))
        plt.colorbar()
        plt.title(f"my domain status: t = {(i * dt):.2f}\n(dt = {dt})", pad=10)
        plt.savefig(f"{wd_my}domains_my{i}.png")
        plt.clf()
        plt.close(plt.gcf())

        plt.imshow(img3, cmap='jet', vmin=np.min(img3), vmax=np.max(img3))
        plt.colorbar()
        plt.title(f"mz domain status: t = {(i * dt):.2f}\n(dt = {dt})", pad=10)
        plt.savefig(f"{wd_mz}domains_mz{i}.png")
        plt.clf()
        plt.close(plt.gcf())

    except FileNotFoundError:
        print(f"File number {i} not found!")
        pass

def parallel_save(max_files , rd):
    Parallel(n_jobs=-1)(delayed(image_save)(i, rd) for i in range(0, max_files + 1, save_after))

def create_animation(images, output_path):
    image_list = [Image.open(img_path).copy() for img_path in images]
    iio.mimsave(output_path, image_list, fps=100)

parallel_save(iteration, read_dir)

print("Images generated!", flush=True)

save_path = Path(wd_mx)
images1 = natsorted(list(save_path.glob('domains_mx*.png')), key=lambda x: x.stem)
save_path = Path(wd_my)
images2 = natsorted(list(save_path.glob('domains_my*.png')), key=lambda x: x.stem)
save_path = Path(wd_mz)
images3 = natsorted(list(save_path.glob('domains_mz*.png')), key=lambda x: x.stem)

create_animation(images1, wd_mx + 'animation.mov')
create_animation(images2, wd_my + 'animation.mov')
create_animation(images3, wd_mz + 'animation.mov')

print("Animation baked!")

stop = time.time()

rt = stop - start
units = 'sec'

if rt > 60 and rt < 3600:
    rt /= 60
    units = 'min'
elif rt > 3600:
    rt /=3600
    units = 'hr'
    
print(f"Script completed in: {rt:.4f} {units}.")
