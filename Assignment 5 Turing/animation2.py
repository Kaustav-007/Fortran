import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
from matplotlib.colors import Normalize

def load_data(filename):
    """Load data for a single species"""
    data = np.loadtxt(filename)
    iterations = data[:, 0].astype(int)
    x_coords = data[:, 1].astype(int)
    y_coords = data[:, 2].astype(int)
    values = data[:, 3]
    
    unique_iters = np.unique(iterations)
    lx, ly = x_coords.max(), y_coords.max()
    
    pattern = np.zeros((len(unique_iters), lx, ly))
    for i, it in enumerate(unique_iters):
        mask = (iterations == it)
        pattern[i, x_coords[mask]-1, y_coords[mask]-1] = values[mask]
    
    return pattern, unique_iters

def save_B_video():
    pattern, iterations = load_data('initializeB.txt')
    
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_title('Turing Pattern - Species B')
    
    # Set up video writer
    writer = FFMpegWriter(fps=1, bitrate=1800)
    
    norm = Normalize(vmin=pattern.min(), vmax=pattern.max())
    im = ax.imshow(pattern[0], cmap='plasma', norm=norm, origin='lower')
    plt.colorbar(im, label='Concentration B')
    
    def update(frame):
        im.set_array(pattern[frame])
        ax.set_title(f'Species B (time {iterations[frame]})')
        return im,
    
    print("Rendering Species B video...")
    with writer.saving(fig, "turing_pattern_B.mp4", dpi=100):
        for frame in range(len(iterations)):
            update(frame)
            writer.grab_frame()
    print("Video saved as 'turing_pattern_B.mp4'")

if __name__ == '__main__':
    save_B_video()
