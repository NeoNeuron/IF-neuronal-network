# Pearson's correlation coefficient are identical to the spike triggered average under my simulation setups
Pearson's correlation coefficient(PCC) measures correlation between spike train and local field potential(LFP) in my program, which can be expressed as:

$$\rho_{xy}(\tau)=\frac{Cov(X(t);Y(t + \tau))}{\sigma_x\sigma_y}$$

where x, representing spike train, is binary sequence, and y represents the LFP. $Cov$ is the notation of covariance and $\sigma$ stands for the standard deviation. For each elements in x, '0' indicates no spiking events happens within time period, dt, which is chosen as a specific value, and otherwise there is a spike. Meanwhile, each elements in y stands for the mean value of continues variable, y, within time period dt.

Expand previous expression as follows:
$$\rho_{xy}(\tau)=\frac{1}{\sigma_x\sigma_y}(E(X(t)Y(t + \tau))-E(X)E(Y))$$
As $\tau$ changes, the second term of the equation doesn't change, which can be threw out as a constant. Meanwhile, X is nonzero only when $t = t_i$, where $t_i$ is the spiking time of i-th spiking events in X. Therefore, PCC can be rewrite as:
$$\rho_{xy}(\tau)=\frac{1}{\sigma_x\sigma_y}\frac{1}{N}\sum_i^nY(t_i + \tau) + C$$
where N is the number of elements in x(and y).

On the other hand, spike triggered average is expressed as:
$$\hat{Y}(\tau) = E(Y(t_i+\tau)) = \frac{1}{n}\sum^n_iY(t_i+\tau)$$
where $t_i$ are the spiking time of i-th spiking events in spike train x, and n are the number of spikes in x.

Comparing the expression of PCC and spike triggered average, we can conclude the that the pattern of curves of them should be identical.
