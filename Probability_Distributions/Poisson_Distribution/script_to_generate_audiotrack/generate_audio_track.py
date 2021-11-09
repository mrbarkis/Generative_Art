import pickle
import matplotlib.pyplot as plt
from scipy.io import wavfile
import numpy as np

with open('collision_times.pickle', 'br') as file:
    t_in_frames = pickle.load(file)

f_s, data = wavfile.read('click.wav')

FPS = 60.0
t_in_seconds = [frame / FPS for frame in t_in_frames]
t_in_samples = [round(t * f_s) for t in t_in_seconds]
max_sample = t_in_samples[-1] + len(data) + 1
print(t_in_frames[-1])
trigger_signal = np.zeros(shape=(max_sample,))
trigger_signal[t_in_samples] = 1.
convolution = np.convolve(data[:,0], trigger_signal, mode='same')

print(len(data))
out = np.concatenate([np.zeros(len(data) // 2), convolution])
out = np.iinfo(np.int16).max * out / np.max(out)

plt.plot(trigger_signal)
plt.plot(out)
plt.show()
print(f"samplerate: {f_s}")
wavfile.write("audio_track.wav", f_s,out.astype(np.int16))
