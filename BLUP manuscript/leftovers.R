\begin{equation} \tag{3.1}\label{eq:3.1}
\begin{gathered}[t]
y_{ij} \sim f \left(\eta^{(y)}_{ij }, \theta^{(y)}_{ij} \right)  \\
g_\eta \left( \eta^{(y)}_{ij} \right) = \mu_0^{(y)} + \mu_j^{(y)}+ \left(\beta_1 + \beta_j^{(y)} \right) x_{ij} \nonumber \\
g_\theta \left( \theta^{(y)}_{ij} \right) = \theta_0^{(y)} + \theta_{j}^{(y)} \nonumber \\
\begin{bmatrix}
\boldsymbol{\alpha}^{(y)} &
\boldsymbol{\beta}^{(y)} &
\boldsymbol{\theta}^{(y)} \end{bmatrix} ^\textrm{T}
 \sim \mathrm{M}\mathrm{VNormal} \left(
\boldsymbol{0},\boldsymbol{\mathrm{P}} \right) \nonumber \\ \nonumber \\
w_{j} \sim \mathrm{Normal} \left( \upmu _{j}, \upsigma _{j} \right) \nonumber \\
\upmu _{j} = \upnu _0 + b _1 \left( \mu^{(y)}_{j} \right) + b _2 \left( \beta_j^{(y)} \right) + 
                b _3 \left( \theta^{(y)}_{j} \right) \\ + \upgamma _4 \left( \mu^{(y)}_j \mu^{(y)}_j \right) + \upgamma _5 \left( \beta^{(y)}_j                             \beta^{(y)}_j \right) + \upgamma _6 \left( \theta^{(y)}_j \theta^{(y)}_j \right) \ \ \ \ \ \  \\ + 
                \upgamma _7 \left( \mu^{(y)}_j \beta^{(y)}_j \right) + \upgamma _8 \left( \mu^{(y)}_j \theta^{(y)}_j \right) + 
                 \upgamma _9 \left( \beta^{(y)}_j \theta^{(y)}_j \right) + \updelta _j \nonumber \\ 
\end{gathered}
\end{equation}





Consider collecting a vector of normally distributed measurements $\mathbf{x}$ for a sample of organisms. For simplicity, we assume that measurements are made with a known measurement error $\sigma$, which could reflect factors such as the calibration error of a hormonal assay, the expected disagreement between trained observers' records, and/or uncertainty due to a prior inferential procedure conducted on the data such as principal component analysis. For any individual $i$ measured in this sample, we know *a priori* that their observed value $x_i$ is not likely to be a reliable, direct measurement of their true trait value. Rather, $x_i$ is just one of many values drawn from a distribution of possible values expected for this animal under an equivalent measurement procedure
$$x_i \sim \textrm{Normal}(\bar{x}_i, \sigma) $$
where $\bar{x}_i$ is the unobserved, true trait value. We'd like to know whether $\mathbf{\bar{x}}$ causes variation in $\mathbf{y}$, which we estimate with the simple regression coefficient
$$\beta_{y \cdot \bar{x}}=\frac{\mathrm{cov}(\mathbf{\bar{x}},\mathbf{y})}{\mathrm{var}(\mathbf{\bar{x}})}$$
We can see that the introduction of measurement error $\mathbf{\epsilon} \sim \textrm{Normal}(0, \sigma)$ into the measurement process, such that  $\mathbf{x}=\mathbf{\bar{x}} + \mathbf{\epsilon}$, $\mathrm{var}(\mathbf{x})$ is expected to increase while $\mathrm{cov}(\mathbf{x},\mathbf{y})$ remains constant. This causes the expected regression coefficient with error-ridden observed values to be attenuated downard toward zero
$$|  {\beta}_{y \cdot x} | < | \beta_{y \cdot \bar{x}} | $$



This is clear enough for a hypothetical case, but one might nonetheless draw the wrong general conclusion from the demonstration. Given that values estimated with error are expected to be downwardly biased when $\beta_{y \cdot \bar{x}} \neq 0$, if one finds empirical evidence that $\hat{\beta}_{y \cdot x} \neq 0$, doesn't this indicate that we should be even *more* confident that a meaningful effect exists? In other words, if we still find support for an association in spite of the random error in our observed measures, isn't that even greater evidence of a robust, true signal beneath the noise? While this is a reasonable first impression, it is important to understand why it is also nonetheless false. Firstly, this conclusion relies on the assumption that the true effect is in fact non-zero, or at least not sufficiently close to zero to be biologically meaningless. Given that empirical research is predicated on the possibility that hypotheses are false and null effects are present, one simply cannot assume that $\hat{\beta}_{y \cdot x} \neq 0$ further indicates that $|  \hat{\beta}_{y \cdot x} | < | \beta_{y \cdot \bar{x}} |$. More, that the use of $dd$, without appropriately accounting for $x$, artificially reduces uncertainty in the estimation of. To see this for a classical estimator, note that ...


Using the observed values $\mathbf{x}$ without accounting for $\sigma$ inhibits robust testing of evolutionary ecological theory- find a clear relationship between $\mathbf{x}$ and some expected outcome $\mathbf {y}$, how can we be confident that this was caused by a failed theoretical prediction rather than random sampling error? 



is that, should a relationship nonetheless be found under these circumstances (e.g. by using a null-hypothesis test), it is an indication that the finding is even more robust as a result. In other words, since 


fallacy Gelman significance under noise.


By ignoring the expected measurement error in a single or average measure, noise would be introduced into any analysis of the association between Thus, in itself, choosing to ignore measurement error in any measurement is a decision to potentially bias statistical inference. 

However, in a perfectly balanced experimental design, it may be relatively straightforward to collect sufficient repeated measurements *r* ${x_{i1},...,x_{ir}}$ so that the true value can be well approximated such that $\frac{1}{r}{}\sum_{i=1}^{r}x_{i} \approx  \bar{x}_i$ with an expected standard error of $\sigma_{\bar{x}_i}=\frac{\sigma}{\sqrt{r}}$. 
things are rarely this simple for behavioral RNs, and as a consequence, the scope for bias is also generally much greater. Firstly, note that we are generally not in a position to know the true value of 


Given these considerations, it's important to address a pontential source

# A novel solution
...

## Selecting model priors
...

# Simulated examples
...