
import math
import pyqtree
import time

class YKJ2TM35:
    def __init__(self):
        self.triangle_points = {}
        self.triangles = {}
        self.triagindexes = {
            'tm35':pyqtree.Index(bbox=(50199.4814, 6582464.0358, 761274.6247, 7799839.8902)),
            'ykj':pyqtree.Index(bbox=(3049997.0487, 6585250.7768, 3761392.2103, 7803169.4895))

        }
        self.loadTriangles()
        self.loadTriangleParams('KKJ_TO_ETRS_TM35FIN.txt','fwd')
        self.loadTriangleParams('ETRS_TM35FIN_TO_KKJ.txt','inv')


    def loadTriangles(self):
        f = open('kkjEUREFFINtriangulationVertices.txt','rb')
        for l in f:
            l = l.strip().split('\t')
            pid = int(l[0])
            self.triangle_points[pid] = {
                'ykj':tuple(reversed(map(float,l[1:3]))),
                'tm35':tuple(reversed(map(float,l[3:5])))
            }
        f.close()

    def loadTriangleParams(self,fn,convtype):
        f = open(fn,'rb')
        for l in f:
            l = l.strip().split('\t')
            l[:4] = map(int,l[:4])
            l[4:] = map(float,l[4:])
            
            tid = l[3]

            if not tid in self.triangles:
                self.triangles[tid] = {
                    'id':tid,
                    'fwd':None,
                    'inv':None,
                    'tids':tuple(l[:3]),
                    'triag': tuple((self.triangle_points[tpid] for tpid in l[:3])),
                    'bbox_ykj':(
                        min((self.triangle_points[tpid]['ykj'][0] for tpid in (l[:3]))),
                        min((self.triangle_points[tpid]['ykj'][1] for tpid in (l[:3]))),
                        max((self.triangle_points[tpid]['ykj'][0] for tpid in (l[:3]))),
                        max((self.triangle_points[tpid]['ykj'][1] for tpid in (l[:3]))),                
                    ),
                    'bbox_tm35':(
                        min((self.triangle_points[tpid]['tm35'][0] for tpid in (l[:3]))),
                        min((self.triangle_points[tpid]['tm35'][1] for tpid in (l[:3]))),
                        max((self.triangle_points[tpid]['tm35'][0] for tpid in (l[:3]))),
                        max((self.triangle_points[tpid]['tm35'][1] for tpid in (l[:3]))),                
                    )


                }

                self.triagindexes['tm35'].insert(item=tid,bbox=self.triangles[tid]['bbox_tm35'])
                self.triagindexes['ykj'].insert(item=tid,bbox=self.triangles[tid]['bbox_ykj'])

            self.triangles[tid][convtype] = {}
            for i,k in enumerate(('a1','a2','dE','b1','b2','dN')):
                self.triangles[tid][convtype][k] = l[4+i]      


        f.close()

    def __doAffineTransformation(self,p,params):
        return (
            params['a1'] * p[1] + params['a2'] * p[0] + params['dE'],
            params['b1'] * p[1] + params['b2'] * p[0] + params['dN']
            )

    def __getTriangleForPoint(self,p,type):

        idxtriags = self.triagindexes[type].intersect((p[0],p[1],p[0],p[1]))

        triag_match = None
        for tid in idxtriags:
            t = self.triangles[tid]
            bbox = t['bbox_%s' % type]
            if p[0] < bbox[0] or p[1] < bbox[1]:
                continue
            if p[0] > bbox[2] or p[1] > bbox[3]:
                continue

            if self.__ptInTriangle(p,(tp[type] for tp in t['triag'])):
                triag_match = t

        return triag_match

    def __ptInTriangle(self,p,triag):
        p0,p1,p2 = triag
        A = 1.0 / 2.0 * (-p1[1] * p2[0] + p0[1] * (-p1[0] + p2[0]) + p0[0] * (p1[1] - p2[1]) + p1[0] * p2[1])
        s = math.copysign((p0[1] * p2[0] - p0[0] * p2[1] + (p2[1] - p0[1]) * p[0] + (p0[0] - p2[0]) * p[1]),A)
        t = math.copysign((p0[0] * p1[1] - p0[1] * p1[0] + (p0[1] - p1[1]) * p[0] + (p1[0] - p0[0]) * p[1]),A)
        return s > 0 and t > 0 and (s + t) < (2.0 * A*math.copysign(1,A))

    def fwd(self,p):
        triag = self.__getTriangleForPoint(p, 'ykj');
        assert triag
        assert triag['fwd']
        return self.__doAffineTransformation(p, triag['fwd']);

if __name__ == '__main__':
    y2e = YKJ2TM35()

