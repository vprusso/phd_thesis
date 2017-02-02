# Ph.D. Thesis : Extended Nonlocal Games
## Version 0.0.1 - 11 November 2016

Software, slides, etc. from my Ph.D. thesis. 

## Abstract

The notions of *entanglement* and *nonlocality* are among the most striking ingredients found in quantum information theory. One tool to better understand these notions is the model of *nonlocal games*; a mathematical framework that abstractly models a physical system. The simplest instance of a nonlocal game involves two players, Alice and Bob, who are not allowed to communicate with each other once the game has started and who play cooperatively against an adversary referred to as the referee. 

The focus of this thesis is a class of games called *extended nonlocal games*, of which nonlocal games are a subset. In an extended nonlocal game, the players initially share a tripartite state *with the referee*. In such games, the winning conditions for Alice and Bob may depend on outcomes of measurements made by the referee, on its part of the shared quantum state, in addition to Alice and Bob's answers to the questions sent by the referee. 

We build up the framework for extended nonlocal games and study their properties and how they relate to nonlocal games. In doing so, we study the types of *strategies* that Alice and Bob may adopt in such a game. For instance, we refer to strategies where Alice and Bob use quantum resources as *standard quantum strategies* and strategies where there is an absence of entanglement as an *unentangled strategy*. These formulations of strategies are purposefully reminiscent of the respective quantum and classical strategies that Alice and Bob use in a nonlocal game, and we also consider other types of strategies with a similar correspondence for the class of extended nonlocal games. 

We consider the *value* of an extended nonlocal game when Alice and Bob apply a particular strategy, again in a similar manner to the class of nonlocal games. Unlike computing the unentangled value where tractable algorithms exist, directly computing the standard quantum value of an extended nonlocal game is an intractable problem. We introduce a technique that allows one to place upper bounds on the standard quantum value of an extended nonlocal game. Our technique is a generalization of what we refer to as the *QC hierarchy* which was studied independently in works by Doherty, Liang, Toner, and Wehner as well as by Navascues, Pironio, and Acin. This technique yields an upper bound approximation for the quantum value of a nonlocal game.

We also consider the question of whether or not the dimensionality of the state that Alice and Bob share as part of their standard quantum strategy makes any difference in how well they can play the game. That is, does there exist an extended nonlocal game where Alice and Bob can win with a higher probability if they share a state where the dimension is infinite? We answer this question in the affirmative and provide a specific example of an extended nonlocal game that exhibits this behavior.   

We study a type of extended nonlocal game referred to as a *monogamy-of-entanglement game*, introduced by Tomamichel, Fehr, Kaniewski, and Wehner, and present a number of new results for this class of game. Specifically, we consider how the standard quantum value and unentangled value of these games relate to each other. We find that for certain classes of monogamy-of-entanglement games, Alice and Bob stand to gain no benefit in using a standard quantum strategy over an unentangled strategy, that is, they perform just as well without making use of entanglement in their strategy. However, we show that there does exist a monogamy-of-entanglement game in which Alice and Bob do perform strictly better if they make use of a standard quantum strategy. We also analyze the *parallel repetition* of monogamy-of-entanglement games; the study of how a game performs when there are multiple instances of the game played independently. We find that certain classes of monogamy-of-entanglement games obey *strong parallel repetition*. In contrast, when Alice and Bob use a non-signaling strategy in a monogamy-of-entanglement game, we find that strong parallel repetition is not obeyed. 


## Dependencies

1. [MATLAB](http://www.mathworks.com/products/matlab/) (free if UW student check IST Software Portal),
2. [CVX](http://cvxr.com/cvx/download/) (free MATLAB package),
3. [QETLAB](http://www.qetlab.com/Main_Page) (free MATLAB package).
4. The software on this page. 

In order to run, place the QETLAB folder, CVX folder, and the content from this repository into a folder in your MATLAB working directory. 
Once you add the folder to your working directory, you are ready to run the examples. 

## Misc. Links:

1. [Nathaniel Johnston QETLAB blog post](http://www.njohnston.ca/2015/04/introducing-qetlab-a-matlab-toolbox-for-quantum-entanglement/) A quick start guide to QETLAB
written by Nathaniel Johnston on his blog. 
2. [How to cite QETLAB](http://www.qetlab.com/How_to_cite)

### References:

[1] "A monogamy of entanglement game with applications to device independent quantum cryptography",
     
      M. Tomamichel, S. Fehr, J. Kaniewski, S. Wehner.,
      
      New Journal of Physics, IOP Publishing, 2013, 15, 103002,
      
      ArXiv: [arxiv:1210.4359](http://arxiv.org/abs/1210.4359)

[2] "Extended nonlocal games and monogamy-of-entanglement games", 

    N. Johnston, R. Mittal, V. Russo, J. Watrous,
    
    Proc. R. Soc. A, 2016, 472, 20160003, 2016,

    ArXiv: [arxiv:1510.02083](https://arxiv.org/abs/1510.02083)
         