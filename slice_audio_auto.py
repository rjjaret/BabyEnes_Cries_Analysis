import os
import librosa
import soundfile as sf
import numpy as np

def detect_cries(y, sr, frame_length=2048, hop_length=512, energy_thresh=0.1, min_silence=0.2, min_cry=0.3):
    energy = librosa.feature.rms(y=y, frame_length=frame_length, hop_length=hop_length)[0]
    frames = np.where(energy > energy_thresh)[0]
    if len(frames) == 0:
        return []
    # Convert frame indices to time
    times = librosa.frames_to_time(frames, sr=sr, hop_length=hop_length)
    # Group contiguous frames
    segments = []
    start = times[0]
    for i in range(1, len(times)):
        if times[i] - times[i-1] > min_silence:
            end = times[i-1]
            if end - start >= min_cry:
                segments.append((start, end, end - start))
            start = times[i]
    # Add last segment
    if times[-1] - start >= min_cry:
        segments.append((start, times[-1], times[-1] - start))
    return segments

def slice_audio_auto(input_path, output_dir, sr=None, **kwargs):
    os.makedirs(output_dir, exist_ok=True)
    y, orig_sr = librosa.load(input_path, sr=sr)
    use_sr = sr or orig_sr
    segments = detect_cries(y, use_sr, **kwargs)
    output_files = []
    base_name = os.path.splitext(os.path.basename(input_path))[0]
    for i, (start, end, _) in enumerate(segments):
        start_sample = int(start * use_sr)
        end_sample = int(end * use_sr)
        cry_audio = y[start_sample:end_sample]
        out_path = os.path.join(output_dir, f"{base_name}_{i+1}.wav")
        sf.write(out_path, cry_audio, use_sr)
        output_files.append(out_path)
    return output_files, segments


def slice_audio_files(input_dir, output_dir, sr=None, **kwargs):
    results = {}
    for root, dirs, files in os.walk(input_dir):
        for filename in files:
            print(f"Processing {filename}...")
            if filename.lower().endswith(('.wav', '.mp3', '.flac')):
                input_path = os.path.join(root, filename)
                # Preserve subdirectory structure in output, but do NOT add filename as a directory
                rel_dir = os.path.relpath(root, input_dir)
                out_subdir = os.path.join(output_dir, rel_dir)
                os.makedirs(out_subdir, exist_ok=True)
                files_out, segs = slice_audio_auto(input_path, out_subdir, sr=sr, **kwargs)
                results[input_path] = (files_out, segs)
    return results

slice_audio_files("data/donateacry_corpus/", "data/donateacry_corpus_sliced/", sr=8000, energy_thresh=0.05)