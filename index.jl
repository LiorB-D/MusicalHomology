using Plots
using WAV
using STFT
using Ripserer
using PersistenceDiagrams


function computePD(fileStr, vol)

    snd, sampFreq = wavread(fileStr)

    N, _ = size(snd)
    rawAmp = snd[1440000:2880000,1] / vol # Only consider 30-minute of song


    W = Int(sampFreq * 0.03)         # Window length
    w = ones(W)       # Rectangular analysis window
    H = 100            # Hop
    L = W - H         # Overlap

    X = stft(rawAmp, w, L)    # Analysis
    s = abs2.(X)         # Compute spectrogram
    println("Fourier Transform complete")

    #heatmap(s) # Display spectrogram


    rawSoundArray = NTuple{4, Float32}[]

    for t=1:length(s[1,:])
        for b=1:length(s[:,1])
            if s[b,t] > 13000
                freq = b * sampFreq / W
                octaveBand = log(2,freq)
                k = floor(octaveBand)
                noteRatio = (freq - 2^k) / (2^(k+1) - 2^k)
                noteX = cos(2pi * noteRatio)
                noteY = sin(2pi * noteRatio)
                v = (t, noteX, noteY, s[b,t] / 1000)
                push!(rawSoundArray, v)
            end
        end
    end
    println("Sound Array Assembled")

    rips = ripserer(rawSoundArray, dim_max = 2, threshold = 3, verbose=true)
    println("Persistence Calculated")
    println(rips)
    for d=1:3
        for h=1:length(rips[d])
            if rips[d][h][2] == Inf
                rips[d][h] = PersistenceInterval(rips[d][h][1], 3)
            end
        end
    end

    return rips
end

april = computePD("April.wav", 1)
blowin = computePD("blowin.wav", 1.5)
youMayBe = computePD("YouMayBe.wav", 1.5)
judas = computePD("Judas.wav", 1)
bleecker = computePD("bleecker.wav", 1)
println("April to blowin")
Bottleneck()(april,blowin)
println("April to You May Be")
Bottleneck()(april,youMayBe)
println("April to Judas") 
Bottleneck()(april, judas) 
println("April to bleecker") 
Bottleneck()(april, bleecker) 
println("blowin to You May Be")
Bottleneck()(blowin,youMayBe)
println("blowin to Judas") 
Bottleneck()(blowin, judas) 
println("blowin to bleecker") 
Bottleneck()(blowin, bleecker) 
println("You May be to Judas") 
Bottleneck()(youMayBe, judas) 
println("You May be to bleecker") 
Bottleneck()(youMayBe, bleecker) 
println("Judas to bleecker") 
Bottleneck()(Judas, bleecker) 

# Frequency of Nth bin - n * sampFreq / window length

