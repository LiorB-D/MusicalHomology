using Plots
using WAV
using STFT
using Ripserer
using PersistenceDiagrams
using Clustering
using Distances



function computePH(fileStr)
    
    snd, sampFreq = wavread(fileStr)

    N, _ = size(snd)
    
    rawAmp = snd[1440000:2880000,1] # Only consider 30-second of song

    # Frequency of Nth bin - n * sampFreq / window length

    W = Int(sampFreq * 0.03)         # Window length
    w = ones(W)       # Rectangular analysis window
    H = 100            # Hop
    L = W - H         # Overlap

    X = stft(rawAmp, w, L)    # Analysis
    s = abs2.(X)         # Compute spectrogram
    println("Fourier Transform complete")
    
    
    #heatmap(s) # Display spectrogram

    
    maxAmp = maximum(s) # Compute maximization

    rawSoundArray = NTuple{4, Float32}[]

    for t=1:length(s[1,:])
        for b=1:length(s[:,1])
            if s[b,t] / maxAmp > 0.025 # Exclude notes that are too insignificant
                freq = b * sampFreq / W
                octaveBand = log(2,freq) #Decimal Octave Band
                k = floor(octaveBand) # Actual Octave Band
                noteRatio = (freq - 2^k) / (2^(k+1) - 2^k) # Position in the Octave Band ie the note
                # Polar Decomposition
                noteX = cos(2pi * noteRatio)
                noteY = sin(2pi * noteRatio)
                normalizedTime = 5 * t / length(s[1,:])
                v = (normalizedTime, noteX, noteY, s[b,t] / 1000)
                push!(rawSoundArray, v)
            end
        end
    end

    println("Sound Array Assembled")
    
    # Cleanup variables
    rawAmp = 0
    w = 0
    X = 0
    s = 0
    # Kmedoids requires dissimilarity matrix
    DissMatrix = pairwise(Euclidean(), rawSoundArray)

    R = kmedoids(DissMatrix, 250)
    DissMatrix = 0 # Cleanup

    println("Medoids Formed")
    rawSoundArray = rawSoundArray[R.medoids] # Reduce the array down to just the medoids
    R = 0
    


    # Running Persistence homology 
    rips = ripserer(rawSoundArray, dim_max = 2, threshold = 1, verbose=true) 
    println("Persistence Calculated")
    println(rips)
    # Ripserer outputs Inf for homology generators that never die, so we need to fix an upper bound to allow the Wasserstein metric computation
    for d=1:3
        for h=1:length(rips[d])
            if rips[d][h][2] == Inf
                rips[d][h] = PersistenceInterval(rips[d][h][1], 3)
            end
        end
    end

    return rips
end



