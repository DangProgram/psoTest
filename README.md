# psoTest

a project of Particle Swarm Optimization

## 一、实验环境

Matlab

## 二、实验原理

粒子群算法（Particle Swarm Optimization,PSO），模拟鸟群觅食行为，每一只鸟（粒子）有当前的位置(Xi)、速度（Vi）以及当前位置的函数适应值（fitness value，即f(xi)），在每一次迭代过程中，每一只鸟更新自己所经历过的个体最优值pbest（即 max{f(xi)}）、更新所有鸟的群体最优值gbest，根据pbest和gbest来决定下次迭代每只鸟的Vi且更新Xi，直到完成迭代。

## 三、实验思考

### （一）迭代过程中参数的影响
<img  src="https://github.com/DangProgram/psoTest/blob/master/pictures/%E5%BE%AE%E4%BF%A1%E6%88%AA%E5%9B%BE_20180921210046.png">

文章中公式（1）（如上图）给出了速度Vi最初的更新公式，最终的到速度大小和方向由三大部分影响，即自身原有的速度、个体最优对该速度、已经群体最优对该速度。可以大体认为（暂时忽视学习因子c1,c2），当个体最优或群体最优离当前位置最远时，速度的改变大；当个体最优和群体最优离当前位置近时，速度的改变小，公式（2）更新当前的位置。
而学习因子c1,c2取值的大小能让速度的更新更加偏向个体的惯性或者个体最优或者群体的最优。选取合适的c1,c2的值，能使得迭代更快，避免陷入局部最优。
但是，公式（1）更新过程中会出现问题，随着不断迭代，Vi是不断增大，但速度过大可能导致无法找到合适的位置（速度太快而错过目标函数的最优）。所以引入惯性因子w（取值为（0,1），一般取0.4）：

![image2](https://github.com/DangProgram/psoTest/blob/master/pictures/%E5%BE%AE%E4%BF%A1%E6%88%AA%E5%9B%BE_20180921212252.png)

引入惯性因子w之后，较好的控制了速度，但是会导致前期搜索速度较慢，因此又引入了w的线性递减权值策略：

![image3](https://github.com/DangProgram/psoTest/blob/master/pictures/%E5%BE%AE%E4%BF%A1%E6%88%AA%E5%9B%BE_20180921212913.png)

当最开始迭代的时候，k=0，wk=wini,当迭代快结束的时候，wk→wend。通过线性递减权值策略，能在迭代过程中更改wk的值，使得开始迭代时，速度保持一个较大值，迭代后期速度逐渐减小，较好的解决了上述问题。
另外，我认为，应当为算法设置结束阈值e，当两次全局最优所对应的值的差的绝对值小于e算法就应该停止，从而减少迭代时间。
###  （二）迭代过程之外参数的影响
迭代之外影响结果的参数主要有粒子个数n和迭代次数，粒子个数或迭代次数太少，都可能无法找到最优值，粒子个数和迭代次数太多又会导致计算资源的浪费和计算时间的增加。
