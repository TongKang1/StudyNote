# Part1. 图论与线性代数

图论基本概念梳理

1.连通性：强连通性（i→j），弱连通性（变成无向图后强连通）和拟强连通（quasi-strongly connected，存在k，使得k→i and k→ j）。显然，拟强连通性与存在一个spanning directed tree（有向生成树）等价。

2.可约：存在排序矩阵（permutation matrix）P，使得 

$$G=P\begin{pmatrix}A&B\\0&C\end{pmatrix}P^T$$

G为原图的邻接矩阵，上述操作的意义为对图的节点（vertex）重新进行编号或排序。显然，不可约与强连通要等价。

3.直径（diameter）：图的直径为任意两个节点之间距离的最大值

$$||(i,j)||_d \leq diameter$$

4.图论中两种重要的矩阵：入射矩阵（interaction matrix）和拉普拉斯矩阵（Laplacian matrix）。A为原图的邻接矩阵，则：

$$A=\begin{pmatrix}0&0.21&3.5\\0&0&2.1\\1.3&0&0\end{pmatrix},  \quad L=\begin{pmatrix}3.71&-0.21&-3.5\\0&2.1&-2.1\\-1.3&0&1.3\end{pmatrix}$$

$$L = D-A, \quad D \; is \; a \; diagonal\;matrix, \;which\;is\;the\;raw\;sum\;of\;A.$$

$$L=MM^T,\;M\;is\;a\;interaction\;matrix.$$

由性质可知，A一定是一个非负矩阵，并且还是对角占优矩阵，按照Gershgorin circle criterion理论，不难看出特征值的分布范围。同时，如果该矩阵不可约，易证：不可约+对角占优矩阵=非奇异矩阵

5.拉普拉斯矩阵的Frobenius normal：

$$L=P\begin{pmatrix}B_1&B_{12}&\cdots&B_{1k}\\&B_2&\cdots&B_{2k}\\ &&\ddots & \vdots \\ &&&B_k\end{pmatrix}P^T,\quad B_i\;are\;square\;irreducible\;matrix.$$

通过该标准型，可以将图的拓扑结构与相应子图之间的连通性揭示出来，可以将一个非强连通的图划分为若干小的强连通的子图去处理。

若把强连通的子图当成一个节点（vertex），就可以得到压缩有向图（condensation directed graph），通过该图，不难找到原图的一个生成有向树（spanning directed tree）。

6.出度与入度：

$$outdegree\;d_o(v),\quad indegree\; d_i(v)$$

$$\delta_o<d_o(v)<\Delta_o,\quad \delta_i<d_i(v)<\Delta_i$$

对于考虑权重的图，这里的度不是简单的输入或输出边的数目，而是输入或输出边的加权和。

7.随机矩阵（stochastic matrices）与双随机矩阵：此处的随机矩阵不是我们通常意义上的随机矩阵，这本书中的随机矩阵指的是行和为1的非负的行随机矩阵（non-negative, row-stochastic matrices）。双随机矩阵则为该矩阵的行与列均满足该条件。

8.随机矩阵的特征值幅角的性质：$\lambda$是该矩阵的特征值，则$1-\lambda$满足

$$\frac{\pi }{n}-\frac{\pi }{2} \leq \phi \leq \frac{\pi }{2}-\frac{\pi }{n}$$

由$A = I - \alpha L$，所以，对于拉普拉斯矩阵，特征值的幅角直接满足上述性质。该结论可以推广到Frobenius标准型下的情况，取决于最大的子块。

$$\frac{\pi }{n_F}-\frac{\pi }{2} \leq \phi \leq \frac{\pi }{2}-\frac{\pi }{n_F}$$

9.m阶可约（m-reducible）: 可以化为如下的Frobenius标准型

$$L=P\begin{pmatrix}B_1&B_{12}&\cdots&&&B_{1k}\\ &\ddots&\cdots&&&B_{2k}\\ &&B_k&B_{k,k+1} & \cdots&B_{k,k+m} \\ &&&B_{k+1}&0&0 \\    &&&&\ddots&0\\ &&&&& B_{k+m} \end{pmatrix}P^T$$

对应到矩阵特征值的情况，L为拉普拉斯矩阵时，特征值0所对应的重数是m。不难看出，将这个矩阵转置（即为逆向图），此时生成有向树的数目至少为m。

所以，1-可约和不可约对应的特征值0的重数都是1，且逆向有向生成树的数目为1个。（why?）

10.左特征值于强连通性：$\omega$为左特征向量，如果$\omega >0,\;s.t.\;\omega^TL=0$，当且仅当强连通时。

代数连通性

$$\lambda_{min}(A) = \lambda_1(A) \leq \lambda_2(A)\leq \cdots \leq \lambda_n(A)=\lambda_{max}(A)$$

1.无向图的连通性：$minimum\;degree\;is\;\delta,\;maximum\;degree\;is\;\Delta$

（1）$\lambda_2\geq0$，不等号在无向图连通时严格成立

（2）$\lambda_2 \leq \frac{n}{n-1}\delta \leq \frac{n}{n-1}\Delta \leq \lambda_n$

（3）if an unweighted graph is not complete, then $\lambda_2 \leq \delta$

（4）$L_K$为完全图的拉普拉斯矩阵，如果$\alpha>0 \; and\; \alpha L_K\;is\;positive\;semidefinite$，则$\lambda_2(L) \geq \frac{n}{\alpha}$

2.有向图的连通性：一共有4个参数用来度量有向图的连通性。

$$K\in \{ x \in \mathbb{R}^V,\;x \perp \textbf{1},\; ||x||=1\}$$

向量$\textbf{1}$为原图的拉普拉斯矩阵特征值为0的特征向量，因此面对强连通图时，K可以剔除为0的那个特征值（存在且唯一）。

（a）$a_1(L)=min_{z\in K}x^TLx=min_{x\neq0,x\perp \textbf{1}} \frac{x^TLx}{x^Tx} =\lambda_2$，这一表达式也被称为瑞利熵。

$b_1(L)=max_{||x||=1} x^TLx = \lambda_{max}$

（1）如果图$\mathcal{G}$的逆向图没有包含一个有向生成树，则$a_1(\mathcal{G}) \leq 0$

（2）如果图不是弱连通，则$a_1(\mathcal{G}) \leq 0$

（3）如果一个图是平衡图：$\mathcal{G}$是强连通图 ↔ $a_1(\mathcal{G}) > 0$

（4）super- and sub-additivity：

$$a_1(\mathcal{G+H}) \geq a_1(\mathcal{G})+ a_1(\mathcal{H})  \\b_1(\mathcal{G+H}) \leq b_1(\mathcal{G})+ b_1(\mathcal{H})$$

（5）

$$a_1(\mathcal{G \times H}) \leq min(a_1(\mathcal{G}), a_1(\mathcal{H})) \leq  max(b_1(\mathcal{G}), b_1(\mathcal{H})) \leq b_1(\mathcal{G \times H})$$

（6）

$$\lambda_1(\frac{1}{2}(L+L^T)) \leq a_1(L) \leq \lambda_2(\frac{1}{2}(L+L^T)), \\ \lambda_{n-1}(\frac{1}{2}(L+L^T)) \leq \lambda_{n}(\frac{1}{2}(L+L^T))$$

（7）如果图$\mathcal{G}$的邻接矩阵A，非对角线上的元素按照$P(A_{ij}=1)=p,\;P(A_{ij}=0)=1-p$，则$a_1(\mathcal{G}) \approx pn ,\;in\;probability\;as\;n\;\to\infty$。这条性质也直观的反映了$a_1(\mathcal{G})$的实际意义。注：这边是依概率收敛！

（b）$a_2(L)=min_{x \in \mathbb{R}^{n-1},\;||Qx||=1} x^TQ^TWLQx= \lambda_{min}(\frac{1}{2}Q^T(WL+L^TW)Q)\\=a_1(WL)$

$W$为一个对角矩阵，满足$\omega^TL=0,\; ||\omega||_{\infty }=1,\;W=diag(\omega)$，且为非负向量；条件$x \in \mathbb{R}^{n-1},\;||Qx||=1$等价于$x \in K$。（显然，W为一个非负的对角矩阵，乘上它以后，不会改变之前特征值的符号，作用类比加权平均？与第一个定义无本质上的不同）

（1）如果图$\mathcal{G}$的逆向图没有包含一个有向生成树，则$a_2(\mathcal{G}) \leq 0$

（2）如果图是强连接，则$0<a_2 \leq Re(\lambda)$，对所有垂直于1向量的特征向量所对应的特征值都成立。

（3）如果图$\mathcal{G}$是强连接的，则

$$a_2(\mathcal{G}) \geq \frac{1-cos(\frac{\pi}{n})} Re(\mathcal{G}) \\ a_2(\mathcal{G}) \geq \frac{C_1e(\mathcal{G})} Re(\mathcal{G})-C_2q \\ where\;r =\frac{max_vw(v)}{min_vw(v)},\; q=max_v{w(v)d_o(v)},\\ C_1=2(cos(\frac{\pi}{n})-cos(\frac{2\pi}{n})), \; C_2=2cos(\frac{\pi}{n})(1-cos(\frac{\pi}{n})),\\ e(\mathcal{G})\;is\;the\;edge\;connectivity.$$

（c）如果图是强连接的，定义：$a_3(L)=min_{x\neq 0,x \perp \textbf{1}} \frac{x^TWLx}{x^T(W- \frac{ww^T}{||w||_1})x}$

如果图是强连接的，则 $a_3 \geq a_2 > 0$

（d）$a_4(L)=min_{1\leq i\leq k} \eta_i,\; where\;\;\eta_i=min_{x\neq 0} \frac{x^TW_iB_ix}{x^TW_ix},\;for\;1\leq i\leq k-1\; \\and\; \eta_k=min_{x\neq 0,x \perp \textbf{1}} \frac{x^TW_kB_kx}{x^T(W_k- \frac{w_kw_k^T}{||w_k||_1})x}$

$a_4>0$当且仅当原图的逆向图含有一个有向生成树。

Q&A：W在一个完整且封闭（不与外界交互）的图内，应该和单位矩阵等价，只有在分析一个非常巨大的图的子图部分时，W才不是单位矩阵，a1与a2才有区别。
