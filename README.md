GUI and main algorithm: discrimination.m

sub-functions for signal quantisation:
    - integration.m: simulate Gaussian noise, signal integration, signal quantisation.
    - signal_peaks.m: events detection, events rating.
        - peaks_processing.m: events rating.

sub-functions for events timing/periodicity:
    - delta_tx.m: compute timing intervals.
    - periodicity.m: linear regression of events timing and computation of normalised MSE. 

sub-functions for minimum variance algorithm:
    - min_variance.m: main script for events discrimination.
    - min_variance_ecg.m: algorithm adapted for ECG signals.
    - min_variance_complexity.m: count the number of operations.

sub-functions for agglomerative clustering/k-means:
    - events_clustering.m: detect events and perform agglomerative clustering.
       - agglo_clustering.m: agglomerative clustering with Ward's criterion.
       - outlier.m: remove clusters of few population. 
    
    - distance.m: city-block metric, euclidian metric.