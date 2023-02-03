using FFTW
using Plots
using WAV
using STFT
using Ripserer


snd, sampFreq = wavread("April.wav")

N, _ = size(snd)
t = 0:1/(N-1):1;
rawAmp = snd[:,1]


W = Int(sampFreq * 0.05)         # Window length
w = ones(W)       # Rectangular analysis window
H = 100            # Hop
L = W - H         # Overlap

X = stft(rawAmp, w, L)    # Analysis
s = abs2.(X)         # Compute spectrogram


#heatmap(s) # Display spectrogram


rawSoundArray = NTuple{5, Float64}[]

for t=1:length(s[1,:])
    for b=1:length(s[:,1])
        if s[b,t] > 500
            freq = b * sampFreq / W
            octaveBand = log(2,freq)
            k = floor(octaveBand)
            noteRatio = (freq - 2^k) / (2^(k+1) - 2^k)
            noteX = cos(2pi * noteRatio)
            noteY = sin(2pi * noteRatio)
            v = (t, octaveBand, noteX, noteY, s[b,t])
            push!(rawSoundArray, v)
        end
    end
end

results = ripserer(rawSoundArray)

plot(results)

# Frequency of Nth bin - n * sampFreq / window length

#=
y = fft(s)

y1 = copy(y)
for i = 1:N
    if abs(y1[i]) > 200
        y1[i] = 0
    end
end

s_new = real(ifft(y1))
wavwrite(s_new, "output1.wav", Fs = sampFreq)

y2 = copy(y)
for i = 1:N
    if abs(y2[i]) <  800
        y2[i] = 0
    end
end

s_new = real(ifft(y2))
wavwrite(s_new, "output2.wav", Fs = sampFreq)

sticks((abs.(y1)))
sticks!((abs.(y2)))

s1,k1 = wavread("output1.wav")
s2,k2 = wavread("output2.wav")

for i = 1:N
    s1[i] += s2[i]
end

wavwrite(s1, "output3.wav", Fs = sampFreq)

=#