using Plots
using WAV
using STFT
using Ripserer
using PersistenceDiagrams
using Clustering
using Distances



function computePD(fileStr)
    
    snd, sampFreq = wavread(fileStr)

    N, _ = size(snd)
    
    rawAmp = snd[1440000:2880000,1] # Only consider 30-second of song


    W = Int(sampFreq * 0.03)         # Window length
    w = ones(W)       # Rectangular analysis window
    H = 100            # Hop
    L = W - H         # Overlap

    X = stft(rawAmp, w, L)    # Analysis
    s = abs2.(X)         # Compute spectrogram
    println("Fourier Transform complete")

    #heatmap(s) # Display spectrogram
    maxAmp = maximum(s)

    rawSoundArray = NTuple{4, Float32}[]

    for t=1:length(s[1,:])
        for b=1:length(s[:,1])
            if s[b,t] / maxAmp > 0.025
                freq = b * sampFreq / W
                octaveBand = log(2,freq)
                k = floor(octaveBand)
                noteRatio = (freq - 2^k) / (2^(k+1) - 2^k)
                noteX = cos(2pi * noteRatio)
                noteY = sin(2pi * noteRatio)
                normalizedTime = 5 * t / length(s[1,:])
                v = (normalizedTime, noteX, noteY, s[b,t] / 1000)
                push!(rawSoundArray, v)
            end
        end
    end
    println("Sound Array Assembled")
    
    rawAmp = 0
    w = 0
    X = 0
    s = 0

    DissMatrix = pairwise(Euclidean(), rawSoundArray)

    R = kmedoids(DissMatrix, 250)
    DissMatrix = 0
    println("Medoids Formed")
    rawSoundArray = rawSoundArray[R.medoids]
    R = 0
    



    rips = ripserer(rawSoundArray, dim_max = 2, threshold = 1, verbose=true)
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

@time april = computePD("April.wav")
@time blowin = computePD("blowin.wav")
@time youMayBe = computePD("YouMayBe.wav")
@time judas = computePD("Judas.wav")
@time bleecker = computePD("bleecker.wav")
println("April to blowin")
Wasserstein()(april,blowin)
println("April to You May Be")
Wasserstein()(april,youMayBe)
println("April to Judas") 
Wasserstein()(april, judas) 
println("April to bleecker") 
Wasserstein()(april, bleecker) 
println("blowin to You May Be")
Wasserstein()(blowin,youMayBe)
println("blowin to Judas") 
Wasserstein()(blowin, judas) 
println("blowin to bleecker") 
Wasserstein()(blowin, bleecker) 
println("You May be to Judas") 
Wasserstein()(youMayBe, judas) 
println("You May be to bleecker") 
Wasserstein()(youMayBe, bleecker) 
println("Judas to bleecker") 
Wasserstein()(judas, bleecker) 

# Frequency of Nth bin - n * sampFreq / window length

