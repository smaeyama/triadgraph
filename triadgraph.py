#!/usr/bin/env python
"""
Draw network of symmetric triad transfer

Functions
---------
(Public use)
triadgraph_symmetric_all
triadgraph_symmetric_kpq
triadgraph_directional_all
triadgraph_directional_kpq
triadgraph_mode2mode_all
symmetrize_triadtransfer
directional_triadtransfer
(Private use)
triadgraph_symmetric_kernel
triadgraph_directional_kernel
triadgraph_mode2mode_kernel
convert_energy2color
"""
import numpy as np
import matplotlib.pyplot as plt
import pygraphviz as pgv
from IPython.display import display, SVG


# Plot symmetric triad transfer S_k^pq
def triadgraph_symmetric_kernel(G,trans,k,p,q,screening,pwidth,nodename):
    """
    Draw network of symmetric triad transfer
    
    Parameters
    ----------
    G : AGraph of pygraphviz
        G = pygraphviz.AGraph(directed=True,strict=False)
    trans : Numpy array
        The symmetric triad transfer function S_k^pq
        Its shape is (n,n,n) where n is the number of modes.
        Its amplitude should be normalized to draw a graph.
    k : int
        index of S_k^pq
    p : int
        index of S_k^pq
    q : int
        index of S_k^pq
    screening : float
        For visibility, draw edges only for |S_k^pq| > screening.
    pwidth : float
        Penwidth for drawing edges, penwidth=pwidth*|S_k^pq|.
    nodename : list
        List of node name, len(nodename) = n where n is the number of modes.
    
    Returns
    -------
    G : AGraph of pygraphviz
        Nodes and edges are added, representing S_k^pq + S_p^qk + S_q^kp = 0
    """
    kpq_junction = "{},{},{}".format(nodename[k],nodename[p],nodename[q])
    wj = np.array([trans[k,p,q], trans[p,q,k], trans[q,k,p]])
    if (k==p and p==q):          # S_k^kk = 0
        pass
    elif (p==q):                 # S_k^pp = - 2*S_p^kp
        if np.abs(wj[0]) > screening: # Screening for visibility
            if wj[0] < 0:
                G.add_edge(nodename[k],kpq_junction,penwidth=pwidth*abs(wj[0]))
                G.add_edge(kpq_junction,nodename[p],penwidth=pwidth*abs(wj[0]))
            else:
                G.add_edge(kpq_junction,nodename[k],penwidth=pwidth*abs(wj[0]))
                G.add_edge(nodename[p],kpq_junction,penwidth=pwidth*abs(wj[0]))
    elif (q==k):                 # S_p^kk + 2*S_k^kp = 0
        if np.abs(wj[1]) > screening: # Screening for visibility
            if wj[1] < 0:
                G.add_edge(nodename[p],kpq_junction,penwidth=pwidth*abs(wj[1]))
                G.add_edge(kpq_junction,nodename[k],penwidth=pwidth*abs(wj[1]))
            else:
                G.add_edge(kpq_junction,nodename[p],penwidth=pwidth*abs(wj[1]))
                G.add_edge(nodename[k],kpq_junction,penwidth=pwidth*abs(wj[1]))
    elif (k==p):                 # S_q^kk + 2*S_k^kq = 0
        if np.abs(wj[2]) > screening: # Screening for visibility
            if wj[2] < 0:
                G.add_edge(nodename[q],kpq_junction,penwidth=pwidth*abs(wj[2]))
                G.add_edge(kpq_junction,nodename[k],penwidth=pwidth*abs(wj[2]))   
            else:
                G.add_edge(kpq_junction,nodename[q],penwidth=pwidth*abs(wj[2]))
                G.add_edge(nodename[k],kpq_junction,penwidth=pwidth*abs(wj[2]))
    else:                        # S_k^pq + S_p^qk + S_q^kp = 0
        if np.abs(2*wj[0]) > screening: # Screening for visibility
            if wj[0] < 0:
                G.add_edge(nodename[k],kpq_junction,penwidth=pwidth*abs(2*wj[0]))
            else:
                G.add_edge(kpq_junction,nodename[k],penwidth=pwidth*abs(2*wj[0]))
        if np.abs(2*wj[1]) > screening: # Screening for visibility
            if wj[1] < 0:
                G.add_edge(nodename[p],kpq_junction,penwidth=pwidth*abs(2*wj[1]))
            else:
                G.add_edge(kpq_junction,nodename[p],penwidth=pwidth*abs(2*wj[1]))
        if np.abs(2*wj[2]) > screening: # Screening for visibility
            if wj[2] < 0:
                G.add_edge(nodename[q],kpq_junction,penwidth=pwidth*abs(2*wj[2]))
            else:
                G.add_edge(kpq_junction,nodename[q],penwidth=pwidth*abs(2*wj[2]))

    return G



def triadgraph_symmetric_all(trans,output=None,title=None,screening=0.1,pwidth=5.0,nodename=None,energy=None):
    """
    Draw network of symmetric triad transfer
    
    Parameters
    ----------
    trans : Numpy array
        The symmetric triad transfer function S_k^pq
        Its shape is (n,n,n) where n is the number of modes.
        Its amplitude should be normalized to draw a graph.
    output : str
        If output == None:
            show a network graph on display.
        else:
            save a png or dot file as path=output.
    title : str, optional
        Title of graph
    screening : float, optional
        For visibility, draw edges only for |S_k^pq| > screening.
        Default: screening=0.1
    pwidth : float, optional
        Penwidth for drawing edges, penwidth=pwidth*|S_k^pq|.
        Default: pwidth=5.0
    nodename : list, optional
        List of node name, len(nodename) = n where n is the number of modes.
    energy : Numpy array, optional
        Energy of the modes
        Its shape is (n) where n is the number of modes.
        Its amplitude should be normalized to draw a graph.
    """
    n = trans.shape[0]
    if title is None:
        G = pgv.AGraph(directed=True,strict=False)
    else:
        G = pgv.AGraph(directed=True,strict=False,label=title)
    if nodename is None:
        nodename = list(range(n))

    # add nodes (Radial layout, color by energy)
    if energy is None:
        energy = np.zeros(n)
    G.add_node(nodename[0],color='red',shape="diamond",pos="0,0",pin=True,\
               style="filled",fillcolor=convert_energy2color(energy[0]))
    for k in range(1,n):
        theta = 2.0*np.pi*(k-1)/(n-1)
        G.add_node(nodename[k],color='red',shape="diamond",pos="{},{}".format(-2*np.sin(theta),2*np.cos(theta)),pin=True,\
                   style="filled",fillcolor=convert_energy2color(energy[k]))

    # add edges
    for k in range(n):
        for p in range(k,n):
            for q in range(p,n):
                triadgraph_symmetric_kernel(G,trans,k,p,q,screening,pwidth,nodename)

    # draw network
    if output is None:
        img = G.draw(prog="fdp", format="svg")#prog=neato|dot|twopi|circo|fdp|nop.
        display(SVG(img))
    elif output[-3:]=="png":
        G.draw(path=output,prog="fdp",format="png")#prog=neato|dot|twopi|circo|fdp|nop.  
    elif output[-3:]=="dot":
        G.draw(path=output,prog="fdp",format="dot")#prog=neato|dot|twopi|circo|fdp|nop.  

    return



def triadgraph_symmetric_kpq(trans,k_in,p_in,q_in,output=None,title=None,screening=0.1,pwidth=5.0,nodename=None,energy=None):
    """
    Draw network of symmetric triad transfer
    
    Parameters
    ----------
    trans : Numpy array
        The symmetric triad transfer function S_k^pq
        Its shape is (n,n,n) where n is the number of modes.
        Its amplitude should be normalized to draw a graph.
    k_in : int
        index of S_k^pq
    p_in : int
        index of S_k^pq
    q_in : int
        index of S_k^pq
    output : str
        If output == None:
            show a network graph on display.
        else:
            save a png or dot file as path=output.
    title : str, optional
        Title of graph
    screening : float, optional
        For visibility, draw edges only for |S_k^pq| > screening.
        Default: screening=0.1
    pwidth : float, optional
        Penwidth for drawing edges, penwidth=pwidth*|S_k^pq|.
        Default: pwidth=5.0
    nodename : list, optional
        List of node name, len(nodename) = n where n is the number of modes.
    energy : Numpy array, optional
        Energy of the modes
        Its shape is (n) where n is the number of modes.
        Its amplitude should be normalized to draw a graph.
    """
    n = trans.shape[0]
    if title is None:
        G = pgv.AGraph(directed=True,strict=False)
    else:
        G = pgv.AGraph(directed=True,strict=False,label=title)
    if nodename is None:
        nodename = list(range(n))
        
    # add nodes (Radial layout, color by energy)
    if energy is None:
        energy = np.zeros(n)
    G.add_node(nodename[0],color='red',shape="diamond",pos="0,0",pin=True,\
               style="filled",fillcolor=convert_energy2color(energy[0]))
    for k in range(1,n):
        theta = 2.0*np.pi*(k-1)/(n-1)
        G.add_node(nodename[k],color='red',shape="diamond",pos="{},{}".format(-2*np.sin(theta),2*np.cos(theta)),pin=True,\
                   style="filled",fillcolor=convert_energy2color(energy[k]))

    # add edges
    triadgraph_symmetric_kernel(G,trans,k_in,p_in,q_in,screening,pwidth,nodename)

    # draw network
    if output is None:
        img = G.draw(prog="fdp", format="svg")#prog=neato|dot|twopi|circo|fdp|nop.
        display(SVG(img))
    elif output[-3:]=="png":
        G.draw(path=output,prog="fdp",format="png")#prog=neato|dot|twopi|circo|fdp|nop.  
    elif output[-3:]=="dot":
        G.draw(path=output,prog="fdp",format="dot")#prog=neato|dot|twopi|circo|fdp|nop.  

    return



# Plot directional representation of symmetric triad transfer D_{k<-q}^p
def triadgraph_directional_kernel(G,trans,k,p,q,screening,pwidth,nodename):
    """
    Draw network of triad transfer based on directional representation
    
    Parameters
    ----------
    G : AGraph of pygraphviz
        G = pygraphviz.AGraph(directed=True,strict=False)
    trans : Numpy array
        Directional representation of the symmetric triad transfer function D_{k<-q}^p
        D_{k<-q}^p is plotted as directional transfer from q to k via tha coupling with p.
        Its shape is (n,n,n) where n is the number of modes.
        Its amplitude should be normalized to draw a graph.
    k : int
        index of D_{k<-q}^p
    p : int
        index of D_{k<-q}^p
    q : int
        index of D_{k<-q}^p
    screening : float
        For visibility, draw edges only for |S_k^pq| > screening.
    pwidth : float
        Penwidth for drawing edges, penwidth=pwidth*|S_k^pq|.
    nodename : list
        List of node name, len(nodename) = n where n is the number of modes.
    
    Returns
    -------
    G : AGraph of pygraphviz
        Edges are added, D_{k<-q}^p = - D_{q<-k}^p
    """
    if (k==p and p==q):
        pass
    elif (p==q):
        wj=trans[k,p,q]
        if np.abs(wj) > screening: # Screening for visibility
            if wj>0:
                G.add_edge(nodename[q],nodename[k],penwidth=pwidth*abs(wj),label=nodename[p])
            else:
                G.add_edge(nodename[k],nodename[q],penwidth=pwidth*abs(wj),label=nodename[p])
    elif (q==k):
        wj=trans[p,q,k]
        if np.abs(wj) > screening: # Screening for visibility
            if wj>0:
                G.add_edge(nodename[k],nodename[p],penwidth=pwidth*abs(wj),label=nodename[q])
            else:
                G.add_edge(nodename[p],nodename[k],penwidth=pwidth*abs(wj),label=nodename[q])
    elif (k==p):
        wj=trans[q,k,p]
        if np.abs(wj) > screening: # Screening for visibility
            if wj>0:
                G.add_edge(nodename[p],nodename[q],penwidth=pwidth*abs(wj),label=nodename[k])
            else:
                G.add_edge(nodename[q],nodename[p],penwidth=pwidth*abs(wj),label=nodename[k])
    else:
        wj=trans[k,p,q]
        if np.abs(wj) > screening: # Screening for visibility
            if wj>0:
                G.add_edge(nodename[q],nodename[k],penwidth=pwidth*abs(wj),label=nodename[p])
            else:
                G.add_edge(nodename[k],nodename[q],penwidth=pwidth*abs(wj),label=nodename[p])
        wj=trans[p,q,k]
        if np.abs(wj) > screening: # Screening for visibility
            if wj>0:
                G.add_edge(nodename[k],nodename[p],penwidth=pwidth*abs(wj),label=nodename[q])
            else:
                G.add_edge(nodename[p],nodename[k],penwidth=pwidth*abs(wj),label=nodename[q])
        wj=trans[q,k,p]
        if np.abs(wj) > screening: # Screening for visibility
            if wj>0:
                G.add_edge(nodename[p],nodename[q],penwidth=pwidth*abs(wj),label=nodename[k])
            else:
                G.add_edge(nodename[q],nodename[p],penwidth=pwidth*abs(wj),label=nodename[k])
    
    return G



def triadgraph_directional_all(trans,output=None,title=None,screening=0.1,pwidth=5.0,nodename=None,energy=None):
    """
    Draw network of triad transfer based on directional representation
    
    Parameters
    ----------
    trans : Numpy array
        Directional representation of the symmetric triad transfer function D_{k<-q}^p
        D_{k<-q}^p is plotted as directional transfer from q to k via tha coupling with p.
        Its shape is (n,n,n) where n is the number of modes.
        Its amplitude should be normalized to draw a graph.
    output : str
        If output == None:
            show a network graph on display.
        else:
            save a png or dot file as path=output.
    title : str, optional
        Title of graph
    screening : float, optional
        For visibility, draw edges only for |S_k^pq| > screening.
        Default: screening=0.1
    pwidth : float, optional
        Penwidth for drawing edges, penwidth=pwidth*|S_k^pq|.
        Default: pwidth=5.0
    nodename : list
        List of node name, len(nodename) = n where n is the number of modes.
    energy : Numpy array, optional
        Energy of the modes
        Its shape is (n) where n is the number of modes.
        Its amplitude should be normalized to draw a graph.
    """
    n = trans.shape[0]
    if title is None:
        G = pgv.AGraph(directed=True,strict=False)
    else:
        G = pgv.AGraph(directed=True,strict=False,label=title)
    if nodename is None:
        nodename = list(range(n))

    # add nodes (Radial layout, color by energy)
    if energy is None:
        energy = np.zeros(n)
    G.add_node(nodename[0],color='red',shape="diamond",pos="0,0",pin=True,\
               style="filled",fillcolor=convert_energy2color(energy[0]))
    for k in range(1,n):
        theta = 2.0*np.pi*(k-1)/(n-1)
        G.add_node(nodename[k],color='red',shape="diamond",pos="{},{}".format(-2*np.sin(theta),2*np.cos(theta)),pin=True,\
                   style="filled",fillcolor=convert_energy2color(energy[k]))

    # add edges
    for k in range(n):
        for p in range(k,n):
            for q in range(p,n):
                triadgraph_directional_kernel(G,trans,k,p,q,screening,pwidth,nodename)

    # draw network
    if output is None:
        img = G.draw(prog="fdp", format="svg")#prog=neato|dot|twopi|circo|fdp|nop.
        display(SVG(img))
    elif output[-3:]=="png":
        G.draw(path=output,prog="fdp",format="png")#prog=neato|dot|twopi|circo|fdp|nop.  
    elif output[-3:]=="dot":
        G.draw(path=output,prog="fdp",format="dot")#prog=neato|dot|twopi|circo|fdp|nop.  
    
    return



def triadgraph_directional_kpq(trans,k_in,p_in,q_in,output=None,title=None,screening=0.1,pwidth=5.0,nodename=None,energy=None):
    """
    Draw network of triad transfer based on directional representation
    
    Parameters
    ----------
    trans : Numpy array
        Directional representation of the symmetric triad transfer function D_{k<-q}^p
        D_{k<-q}^p is plotted as directional transfer from q to k via tha coupling with p.
        Its shape is (n,n,n) where n is the number of modes.
        Its amplitude should be normalized to draw a graph.
    k_in : int
        index of D_{k<-q}^p
    p_in : int
        index of D_{k<-q}^p
    q_in : int
        index of D_{k<-q}^p
    output : str
        If output == None:
            show a network graph on display.
        else:
            save a png or dot file as path=output.
    title : str, optional
        Title of graph
    screening : float, optional
        For visibility, draw edges only for |S_k^pq| > screening.
        Default: screening=0.1
    pwidth : float, optional
        Penwidth for drawing edges, penwidth=pwidth*|S_k^pq|.
        Default: pwidth=5.0
    nodename : list
        List of node name, len(nodename) = n where n is the number of modes.
    energy : Numpy array, optional
        Energy of the modes
        Its shape is (n) where n is the number of modes.
        Its amplitude should be normalized to draw a graph.
    """
    n = trans.shape[0]
    if title is None:
        G = pgv.AGraph(directed=True,strict=False)
    else:
        G = pgv.AGraph(directed=True,strict=False,label=title)
    if nodename is None:
        nodename = list(range(n))

    # add nodes (Radial layout, color by energy)
    if energy is None:
        energy = np.zeros(n)
    G.add_node(nodename[0],color='red',shape="diamond",pos="0,0",pin=True,\
               style="filled",fillcolor=convert_energy2color(energy[0]))
    for k in range(1,n):
        theta = 2.0*np.pi*(k-1)/(n-1)
        G.add_node(nodename[k],color='red',shape="diamond",pos="{},{}".format(-2*np.sin(theta),2*np.cos(theta)),pin=True,\
                   style="filled",fillcolor=convert_energy2color(energy[k]))

    # add edges
    triadgraph_directional_kernel(G,trans,k_in,p_in,q_in,screening,pwidth,nodename)

    # draw network
    if output is None:
        img = G.draw(prog="fdp", format="svg")#prog=neato|dot|twopi|circo|fdp|nop.
        display(SVG(img))
    elif output[-3:]=="png":
        G.draw(path=output,prog="fdp",format="png")#prog=neato|dot|twopi|circo|fdp|nop.  
    elif output[-3:]=="dot":
        G.draw(path=output,prog="fdp",format="dot")#prog=neato|dot|twopi|circo|fdp|nop.
    
    return



# Plot the contracted mode-to-mode transfer D_{k<-q}
def triadgraph_mode2mode_kernel(G,mode2mode,k,q,screening,pwidth,nodename):
    """
    Draw network of mode-to-mode transfer contracted over the index of mediator
      D_{k<-q} = \sum_p D_{k<-q}^p
    
    Parameters
    ----------
    G : AGraph of pygraphviz
        G = pygraphviz.AGraph(directed=True,strict=False)
    mode2mode : Numpy array
        Mode-to-mode transfer contracted over the index of mediator D_{k<-q}.
        Its shape is (n,n) where n is the number of modes.
        Its amplitude should be normalized to draw a graph.
    k : int
        index of D_{k<-q}
    q : int
        index of D_{k<-q}
    screening : float
        For visibility, draw edges only for |D_{k<-q}| > screening.
    pwidth : float
        Penwidth for drawing edges, penwidth=pwidth*|D_{k<-q}|.
    nodename : list
        List of node name, len(nodename) = n where n is the number of modes.
    
    Returns
    -------
    G : AGraph of pygraphviz
        Edges are added, D_{k<-q} = - D_{q<-k}
    """
    if (k==q):
        pass
    else:
        wj=mode2mode[k,q]
        if np.abs(wj) > screening: # Screening for visibility
            if wj>0:
                G.add_edge(nodename[q],nodename[k],penwidth=pwidth*abs(wj))
            else:
                G.add_edge(nodename[k],nodename[q],penwidth=pwidth*abs(wj))
    
    return G



def triadgraph_mode2mode_all(mode2mode,output=None,title=None,screening=0.1,pwidth=5.0,nodename=None,energy=None):
    """
    Draw network of mode-to-mode transfer contracted over the index of mediator
      D_{k<-q} = \sum_p D_{k<-q}^p
    
    Parameters
    ----------
    mode2mode : Numpy array
        Mode-to-mode transfer contracted over the index of mediator D_{k<-q}.
        Its shape is (n,n) where n is the number of modes.
        Its amplitude should be normalized to draw a graph.
    output : str
        If output == None:
            show a network graph on display.
        else:
            save a png or dot file as path=output.
    title : str, optional
        Title of graph
    screening : float, optional
        For visibility, draw edges only for |D_{k<-q}| > screening.
        Default: screening=0.1
    pwidth : float, optional
        Penwidth for drawing edges, penwidth=pwidth*|D_{k<-q}|.
        Default: pwidth=5.0
    nodename : list
        List of node name, len(nodename) = n where n is the number of modes.
    energy : Numpy array, optional
        Energy of the modes
        Its shape is (n) where n is the number of modes.
        Its amplitude should be normalized to draw a graph.
    """
    n = mode2mode.shape[0]
    if title is None:
        G = pgv.AGraph(directed=True,strict=False)
    else:
        G = pgv.AGraph(directed=True,strict=False,label=title)
    if nodename is None:
        nodename = list(range(n))

    # add nodes (Radial layout, color by energy)
    if energy is None:
        energy = np.zeros(n)
    G.add_node(nodename[0],color='red',shape="diamond",pos="0,0",pin=True,\
               style="filled",fillcolor=convert_energy2color(energy[0]))
    for k in range(1,n):
        theta = 2.0*np.pi*(k-1)/(n-1)
        G.add_node(nodename[k],color='red',shape="diamond",pos="{},{}".format(-2*np.sin(theta),2*np.cos(theta)),pin=True,\
                   style="filled",fillcolor=convert_energy2color(energy[k]))

    # add edges
    for k in range(n):
        for q in range(k,n):
            triadgraph_mode2mode_kernel(G,mode2mode,k,q,screening,pwidth,nodename)

    # draw network
    if output is None:
        img = G.draw(prog="fdp", format="svg")#prog=neato|dot|twopi|circo|fdp|nop.
        display(SVG(img))
    elif output[-3:]=="png":
        G.draw(path=output,prog="fdp",format="png")#prog=neato|dot|twopi|circo|fdp|nop.  
    elif output[-3:]=="dot":
        G.draw(path=output,prog="fdp",format="dot")#prog=neato|dot|twopi|circo|fdp|nop.  
    
    return



def convert_energy2color(energy):
    """
    Coloring node by its energy
    
    Parameter
    ---------
    energy : Numpy array
        Energy of the node
    
    Returns
    -------
    fillcolor : str
        Color of the node in RBGA
    """
    if energy==0:
        fillcolor="#00000000"
    else:
        cm = plt.get_cmap("OrRd",256)
        r,g,b,a=cm(min(int(255*energy),255))
        #r,g,b,a=int(255*r),int(255*g),int(255*b),int(255*a)
        r,g,b,a=int(255*r),int(255*g),int(255*b),int(255*a*energy)
        fillcolor="#{:02x}{:02x}{:02x}{:02x}".format(r,g,b,a)
    return fillcolor



# Calculate symmetric triad transfer S_k^pq
def symmetrize_triadtransfer(trans,time_axis=3):
    """
    Symmetrize triad transfer, S_k^pq = S_k^qp
    
    From assymetric transfer A_k^pq,
        S_k^pq = 0.5 * (A_k^pq + A_k^qp)
    
    Parameters
    ----------
    trans : Numpy array
        Assymetric triad transfer function A_k^pq
        Its shape is,
            (ntime,nk,np,nq) when time_axis=0
            (nk,np,nq,ntime) when time_axis=3
        where nk=np=nq is the number of modes, and ntime is the number of time steps.
    time_axis : int, optional
        The axis along temporal variation of trans
        time_axis should be 0 or 3. Default time_axis=3.
            
    Returns
    -------
    symmetric_trans : Numpy array
        Symmetric triad transfer function S_k^pq
        
    Theory
    ------
    * Symmetry
        S_k^pq = S_k^qp
    * Detailed balance (conservation law among triad)
        S_k^pq + S_p^qk + S_q^kp = 0
    """
    if (time_axis==0):
        symmetric_trans = 0.5 * (trans + np.transpose(trans,axes=(0,1,3,2))) # axes=(2,3) are (p,q)
    elif (time_axis==3):
        symmetric_trans = 0.5 * (trans + np.transpose(trans,axes=(0,2,1,3))) # axes=(1,2) are (p,q)
    else:
        print("# time_axis should be 0 or 3. time_axis = ", time_axis)

    return symmetric_trans



# Calculate directional representation D_{k<-q}^p (from symmetric triad transfer S_k^pq) 
def directional_triadtransfer(symmetric_trans,time_axis=3):
    """
    Directional representation of symmetric triad transfer, D_{k<-q}^p
    D_{k<-q}^p represents symmetric triad transfer as directional transfer from q to k via tha coupling with p.
    
    
    Parameters
    ----------
    symmetric_trans : Numpy array
        Symetric triad transfer function S_k^pq
        Its shape is,
            (ntime,nk,np,nq) when time_axis=0
            (nk,np,nq,ntime) when time_axis=3
        where nk=np=nq is the number of modes, and ntime is the number of time steps.
    time_axis : int, optional
        The axis along temporal variation of trans
        time_axis should be 0 or 3. Default time_axis=3.
            
    Returns
    -------
    directional_trans : Numpy array
        Directional representation D_{k<-q}^p
    
    Theory
    ------
    Because of detailed balance of symmetric triad transfer, S_k^pq + S_p^qk + S_q^kp = 0,
    the interactions among (k,p,q) are always 2 giver(+) - 1 taker(-) or 1 giver(+) - 2 taker(-).
    
    The directional representation splits the symmetric triad transfer according to 
    the rule of "No simultanous gain/loss", or equivalently, 
    "Giver should give, taker should take", or "Minimize |D_{k<-q}^p|+|D_p^qk|+|D_q^kp|".
    
    * No simulaneous gain/loss rule
        If sign(S_k^pq) == sign(S_q^kp), which means both k and p are givers (or takers), then set
        D_{k<-q}^p = 0

    * Conservation law between two (k,q) via a mediator (p)
        D_{k<-q}^p = - D_{q<-k}^p
        
    * Relation with symmetric triad transfer
        S_k^pq = 0.5 * (D_{k<-q}^p + D_{k<-p}^q)
        S_p^qk = 0.5 * (D_{p<-k}^q + D_{p<-q}^k)
        S_q^kp = 0.5 * (D_{q<-p}^k + D_{q<-k}^p)
    """
    n=symmetric_trans.shape[1]
    if (time_axis==0):
        ntime=symmetric_trans.shape[0]
    elif (time_axis==3):
        ntime=symmetric_trans.shape[3]
    else:
        print("# time_axis should be 0 or 3. time_axis = ", time_axis)

    directional_trans=np.zeros_like(symmetric_trans)
    wditr=np.zeros((n,n,n))
    for t in range(ntime):
        if (time_axis==0):
            wtr=symmetric_trans[t,:,:,:]
        elif (time_axis==3):
            wtr=symmetric_trans[:,:,:,t]
        for k in range(n):
            for p in range(k,n):
                for q in range(p,n):
                    wj = np.array([wtr[k,p,q], wtr[p,q,k], wtr[q,k,p]])
                    if (k==p and p==q):          # S_k^kk = 0
                        wditr[k,p,q] = 0.0
                    elif (k==p or p==q or q==k): # S_k^pp = - 2*S_p^kp
                        if (np.prod(wj) > 0):    ## 2 giver(-) 1 taker(+)
                            arg=np.argmax(wj)
                            if arg==0:           ### k is the taker, p=q is giver
                                wt=k
                                wg=p
                            elif arg==1:         ### p is the taker, k=q is giver
                                wt=p
                                wg=q
                            else:                ### q is the taker, k=p is giver
                                wt=q
                                wg=k
                            wditr[wg,wt,wg] = 0.0           # Giver should only give.
                            wditr[wg,wg,wt] = 2*wtr[wg,wg,wt] # giver -> taker (positive value)
                            wditr[wt,wg,wg] = - wditr[wg,wg,wt]
                        else:                    ## 1 giver(-) 2 taker(+)
                            arg=np.argmin(wj)
                            if arg==0:           ### k is the giver, p=q is taker
                                wg=k
                                wt=p
                            elif arg==1:         ### p is the giver, k=q is taker
                                wg=p
                                wt=q
                            else:                ### q is the giver, k=p is taker
                                wg=q
                                wt=k
                            wditr[wt,wg,wt] = 0.0           # Taker should only take.
                            wditr[wt,wt,wg] = 2*wtr[wt,wt,wg] # taker -> giver (negative value)
                            wditr[wg,wt,wt] = - wditr[wt,wt,wg]
                    else:                        # S_k^pq + S_p^qk + S_q^kp = 0
                        if (np.prod(wj) > 0):    ## 2 giver(-) 1 taker(+)
                            arg=np.argmax(wj)
                            if arg==0:           ### k is the taker, p,q are givers
                                wt=k
                                wg1=p
                                wg2=q
                            elif arg==1:         ### p is the taker, q,k are givers
                                wt=p
                                wg1=q
                                wg2=k
                            else:                ### q is the taker, k,p are givers
                                wt=q
                                wg1=k
                                wg2=p
                            wditr[wg1,wt,wg2] = 0.0               # Givers should only give.
                            wditr[wg2,wt,wg1] = 0.0               # Givers should only give.
                            wditr[wg1,wg2,wt] = 2*wtr[wg1,wg2,wt] # taker -> giver1 (negative value)
                            wditr[wg2,wg1,wt] = 2*wtr[wg2,wg1,wt] # taker -> giver2 (negative value)
                            wditr[wt,wg2,wg1] = - wditr[wg1,wg2,wt]
                            wditr[wt,wg1,wg2] = - wditr[wg2,wg1,wt]
                        else:                    ## 1 giver(-) 2 taker(+)
                            arg=np.argmin(wj)
                            if arg==0:           ### k is the giver, p,q are takers
                                wg=k
                                wt1=p
                                wt2=q
                            elif arg==1:         ### p is the giver, q,k are takers
                                wg=p
                                wt1=q
                                wt2=k
                            else:                ### q is the giver, k,p are takers
                                wg=q
                                wt1=k
                                wt2=p
                            wditr[wt1,wg,wt2] = 0.0               # Takers should only take.
                            wditr[wt2,wg,wt1] = 0.0               # Takers should only take.
                            wditr[wt1,wt2,wg] = 2*wtr[wt1,wt2,wg] # giver -> taker1 (positive value)
                            wditr[wt2,wt1,wg] = 2*wtr[wt2,wt1,wg] # giver -> taker2 (positive value)
                            wditr[wg,wt2,wt1] = - wditr[wt1,wt2,wg]
                            wditr[wg,wt1,wt2] = - wditr[wt2,wt1,wg]
        if (time_axis==0):
            directional_trans[t,:,:,:] = wditr
        elif (time_axis==3):
            directional_trans[:,:,:,t] = wditr
    
    return directional_trans
