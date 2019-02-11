# Gaussian Graphical model
Model:
\begin{align}
\begin{split}
{\bf Y}_{t} = \rho {\bf Y}_{t-1} + \bm{\epsilon}_t, t=2,3,\ldots,T\\
\bm{\epsilon}_t \sim {\mathcal{N}}_p\left(0,(LD^{-1}L^T)^{-1}\right)\\
L_{ij}|D_{jj} \sim q\mathbbm{1}_{\{0\}}+(1-q) \mathcal{N}\left(0,D_{jj}\tau^2\right)\\
D_{jj} \sim \text{Inverse-Gamma}(\alpha,\beta)\\ 
\rho \sim 1 \text(improper prior)
\end{split}
\end{align}

Let $S=\sum\limits_{t=2}^{T}({\bf Y}_t-\rho {\bf Y}_{t-1})({\bf Y}_t-\rho {\bf} Y_{t-1})^T$. The joint density of the parameter vector $(L,D,\rho)$ conditioned on data ${\bf Y}$ is given by
\begin{align}\label{eq1:joint}
\pi(L,D,\rho|{\bf Y}) \propto &\left(\prod\limits_{j=1}^p D_{jj}^{-\frac{n-1}2}\exp\left(-\frac {L_{.j}^TSL_{.j} }{2  D_{jj}}\right)\right) \nonumber \\
& \times \prod\limits_{j=1}^{p-1}\prod\limits_{i=j+1}^{p}\left[q \mathbbm{1}_{\{0\}}+(1-q)(2\pi D_{jj}\tau^{2})^{-\frac 1 2}\exp\left(-\frac{L_{ij}^{2}}{2\tau^2D_{jj}}\right)\right] \prod\limits_{j=1}^p D_{jj}^{-\alpha-1}\exp\left(-D_{jj}^{-1}\beta\right)
\end{align}
Based on \ref{eq1:joint}, the following conditional distributions of $L_{ij}$ and $D_{jj}$ can be derived in straightforward fashion.
\begin{align}
L{ij} ~|~ L_{-ij},D_{jj},\rho,{\bf Y} \sim q^* \mathbbm{1}_{\{0\}} +(1-q^*) \mathcal{N}\left(-\frac{b_{ij}}{a_{ij}},\frac{D_{jj}}{a_{ij}}\right),
\end{align}
where $a_{ij}=S_{ii}+\frac{1}{\tau^2}$, $b_{ij}=\sum\limits_{j\leq l \leq p, l\neq i} L_{lj}S_{li}$, $ q^*=\left(1+\frac{1-q}{q}\left(\tau^2a_{ij}\right)^{-\frac 1 2}\exp\left\{\frac {b_{ij}^2} {2 D_{jj}a_{ij}}\right\}\right)^{-1}$.
\begin{align}
D_{jj} ~|~ L_{.j}, \rho, {\bf Y} \sim \mathrm{Inverse-Gamma}(\alpha_{j},\beta_j),
\end{align}
where $\alpha_j = \frac{2\alpha+p-j-k+n-1}{2}$, $\beta_j=\frac 1 2\left( L_{.j}^TSL_{.j}+\sum\limits_{i=j+1}^p \frac{L_{ij}^2}{\tau^2}+2\beta\right)$, $k = \abs{\mathcal{L}_j^{0}=\left\{L_{ij}=0, j+1\leq i \leq p \right\}}$.
