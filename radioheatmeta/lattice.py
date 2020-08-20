#  Copyright (C) 2020 Yang Long (longyang_123@yeah.net)
# 
#  RadioHeatMeta is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  RadioHeatMeta is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

import numpy as np 

DBL_EPSILON = 2.2204460492503131e-16 # Smallest double value

class Lattice:
    '''
    Lattice for multilayer system
    '''
    def __init__(self, ax=[1,0], ay=[0,1], dimension="no", verbose=False):
        self.__lattice = [np.array(ax), np.array(ay)] 

        self.__angle = 0
        self.__area = 0
        self.__reciprocal_lattice = None

        self.__dimension = dimension

        self.verbose = verbose

        self.__Gx = None
        self.__Gy = None
        self.__nG = None

        self.__init_lattice()

    def __init_lattice(self):
        '''
        Initiate the key parameters of the lattice system
        '''
        ax, ay = self.__lattice
        self.__area = np.abs(ax[0]*ay[1] - ax[1]*ay[0]) # TODO may bug for 1D case
        self.__angle = np.arccos(
            np.dot(ax, ay) / (np.square(np.dot(ax, ax)) * np.square(np.dot(ay, ay)))
        )

        # R: the 90 degree rotation matrix
        R = np.array([
            [0 , -1],
            [1 , 0]
        ])

        bx = 2*np.pi* R.dot(ay) / (ax.dot(R.dot(ay)))
        by = 2*np.pi* R.dot(ax) / (ay.dot(R.dot(ax)))
        self.__reciprocal_lattice = [bx, by]

    def get_lattice(self):
        return self.__lattice

    def get_dimension(self):
        return self.__dimension

    def get_reciprocal_lattice(self):
        return self.__reciprocal_lattice

    def get_angle(self):
        return self.__angle

    def get_area(self):
        return self.__area

    def get_nG(self):
        if (self.__nG is None):
            raise Exception("Please initiate the reciprocal G matrix firstly!")
        return self.__nG

    def get_G(self):
        if (self.__Gx is None):
            raise Exception("Please initiate the reciprocal G matrix firstly!")
        return self.__Gx, self.__Gy

    def init_G(self, nG, truncation="parallelogramic"):
        '''
        Compute G matrix 
        '''
        if self.verbose:
            print("Initiating G matrix with the dimension [{dimension}] and the truncation [{truncation}] ... ".format(
                dimension=dimension,
                truncation=truncation
            ))

        dimension = self.__dimension

        if (nG <=0) and (dimension != "no"):
            raise Exception("Require number of G more than 1!")
        if dimension not in ["no", "one", "two"]:
            raise Exception("Unknown dimension type: {dimension}. Only support: no, one or two.".format(dimension=dimension))
        if truncation not in ["circular", "parallelogramic"]:
            raise Exception("Unknown truncation type: {truncation}. Only support: circular, parallelogramic".format(truncation=truncation))

        if (dimension == "no"):
            Gx, Gy, nG = __get_G_parallelogramic(nG=1, dimension)
        elif (dimension == "one"):
            Gx, Gy, nG = __get_G_parallelogramic(nG, dimension)
        else: # 2D dimension
            if (truncation == "circular"):
                Gx, Gy, nG = __get_G_circular(nG, dimension)
            else:
                Gx, Gy, nG = __get_G_parallelogramic(nG, dimension)

        if self.verbose:
            print("G matrix with the number [{nG}]".format(
                nG=nG
            ))

        self.__Gx = Gx
        self.__Gy = Gy
        self.__nG = nG

    def __is_G_same(self, Gi, Gj, Lk):
        Lk_prod = [
            Lk[0]*Lk[0] + Lk[1]*Lk[1],
            2*(Lk[0]*Lk[2] + Lk[1]*Lk[3]),
            Lk[2]*Lk[2] + Lk[3]*Lk[3]
        ]
        ilen = np.hypot(
            Gi[0]*Lk[0] + Gi[1]*Lk[2],
            Gi[0]*Lk[1] + Gi[1]*Lk[3]        
        )
        jlen = np.hypot(
            Gj[0]*Lk[0] + Gj[1]*Lk[2],
            Gj[0]*Lk[1] + Gj[1]*Lk[3]
        )
        maxlen = max(ilen, jlen)

        dv = [
            Gi[0]*Gi[0] - Gj[0]*Gj[0],
            Gi[0]*Gi[1] - Gj[0]*Gj[1],
            Gi[1]*Gi[1] - Gj[1]*Gj[1]
        ]
        d = abs(dv[0]*Lk_prod[0] + dv[1]*Lk_prod[1] + dv[2]*Lk_prod[2])

        return d < 2*DBL_EPSILON*maxlen

    def __sort_G(self, G, Lk_prod):
        '''
        Sort G based on the length
        '''
        u2, uv2, v2 = Lk_prod
        Glen = np.zeros(len(G))
        for i in range(len(G)):
            Gi = G[i]
            Glen[i] = np.sqrt(Gi[0]*Gi[0]*u2 + Gi[0]*Gi[1]*uv2 + Gi[1]*Gi[1]*v2)
        
        sorted_indexs = np.argsort(Glen)
        Gnew = G[sorted_indexs]
        return Gnew

    def __get_G_circular(self, nG, dimension):
        '''
        Compute G matrix for the system using circular truncation
        '''
        bx, by = self.__reciprocal_lattice
        Lk = [
            bx[0]/(2*np.pi),
            bx[1]/(2*np.pi),
            by[0]/(2*np.pi),
            by[1]/(2*np.pi)
        ]

        u = np.hypot(Lk[0], Lk[1])
        v = np.hypot(Lk[2], Lk[3])
        u2 = np.square(u)
        v2 = np.square(v)
        uv = Lk[0]*Lk[2] + Lk[1]*Lk[3]
        Lk_prod = [u2, 2*uv, v2]
        uxv = np.abs(Lk[0]*Lk[3] - Lk[1]*Lk[2])
        circ_area = nG*uxv
        circ_radius = np.sqrt(circ_area/np.pi) + u + v

        u_extent = 1 + int(circ_radius/(u*np.sqrt(1-uv*uv/(u2*v2))))
        v_extent = 1 + int(circ_radius/(v*np.sqrt(1-uv*uv/(u2*v2))))
        uext21 = 2*u_extent + 1
        vext21 = 2*v_extent + 1
        G_temp = np.zeros((uext21*vext21,2))

        for i in range(uext21):
            for j in range(vext21):
                G_temp[i+j*uext21, 0] = i - u_extent
                G_temp[i+j*uext21, 1] = j - v_extent

        G_temp = self.__sort_G(G_temp, Lk_prod)

        length = uext21*vext21
        if (nG < length):
            length = nG
        for i in range(length,0,-1):
            if not self.__is_G_same(G_temp[i-1], G_temp[i], Lk):
                # Get the circular radius in the momentum space
                nG = i
    
        Gx = G_temp[:nG, 0]
        Gy = G_temp[:nG, 1]

        bx, by = self.__reciprocal_lattice
        Gx_temp = Gx*bx[0] + Gy*by[0]
        Gy = Gx*bx[1] + Gy*by[1]
        Gx = Gx_temp
        return Gx, Gy, nG

    def __get_G_parallelogramic(self, nG, dimension):
        '''
        Compute G matrix for the system using parallelogramic trunctions
        '''
        if (dimension == "one"): # 1D case
            m = int(nG/2)
            nG = 2*m + 1
            Gx_list, Gy_list = np.zeros(nG), np.zeros(1)
            for i in range(-M, M+1):
                Gx_list(i+M) = i
            Gx_mat, Gy_mat = np.meshgrid(Gx_list, Gy_list)
        else: # 2D case
            n_root = int(np.sqrt(nG))
            if (n_root % 2 ==0) and (n_root > 0):
                n_root -= 1
            m = int(n_root / 2)

            G_list = np.zeros(n_root)
            for i in range(-M, M+1):
                G_list(i+M) = i
            Gx_mat, Gy_mat = np.meshgrid(G_list, G_list)
            nG = np.power(n_root, 2)

        bx, by = self.__reciprocal_lattice
        Gx_mat_temp = bx[0] * Gx_mat + by[0] * Gy_mat
        Gy_mat = bx[1] * Gx_mat + by[1] * Gy_mat
        Gx_mat = Gx_mat_temp
        return Gx_mat, Gy_mat, nG