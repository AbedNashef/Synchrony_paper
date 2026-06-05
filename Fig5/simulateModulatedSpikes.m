function spikeTrains = simulateModulatedSpikes(realTrace, fs, baselineRate, modulationStrength, nNeurons, dt)
% realTrace: vector of your real trace (e.g., stimulus), size [1 x T]
% fs: sampling rate of realTrace (Hz)
% baselineRate: constant baseline firing rate (Hz)
% modulationStrength: scalar to tune modulation effect
% nNeurons: number of neurons to simulate
% dt: time step (s), e.g., 1/fs

T = length(realTrace);
t = (0:T-1) * dt;

% Normalize and scale the real trace to a modulation in rate
modSignal = (realTrace - mean(realTrace)) / std(realTrace);
% instantaneous rate across time (Hz)
instantRate = baselineRate * (1 + modulationStrength * modSignal);
instantRate(instantRate < 0) = 0; % no negative rates

spikeTrains = zeros(nNeurons,length(realTrace));

for neuronIdx = 1:nNeurons
    % For each time bin, draw from Poisson
    lambda = instantRate * dt; % probability per dt
    spikes = rand(1, T) < lambda;
    spikeTrains(neuronIdx,:) = spikes;  % spike times (s)
end
end