#ifndef FEMOBJECT_H
#define FEMOBJECT_H


#include <trianglesystem/gmtrianglesystem>
#include<math.h>
#include<core/gmarraylx>


using namespace std;

class FemObj: public GMlib::TriangleFacets<float> {

GM_SCENEOBJECT(FemObj)
public:
   FemObj(){  //initialization
    _force=0;
    _swp=true;
}


   ~FemObj(){}


//-------- random points in drum

   void randomTriangulate(float radius, int triangles)
   {
    double s=std::sqrt(( 4*M_PI*radius*radius)/(triangles*std::sqrt(3)));// triangle edge's lengths
    //cout<<"length of s"<<s<<endl;
    int bpoints =( M_2PI*radius)/s ;
    //cout<<"number of boundary:"<<bpoints<<endl;
    int points=(triangles+2+bpoints)/2;
    //cout<<"number of inner points:"<<points<<endl;

    auto _vec=GMlib::Vector<float,2>(radius,0);

    for (int i=0; i<bpoints; i++ )
    {

    GMlib::SqMatrix<float,2> rotationMatrix(GMlib::Angle(M_2PI/bpoints), GMlib::Vector<float,2>(1.0f,0.0f),GMlib::Vector<float,2>(0.0f,1.0f));

    _vec=(rotationMatrix *_vec);
    //cout<<_vec<<endl;
    this->insertAlways(GMlib::TSVertex<float>(static_cast<GMlib::Point<float,2>>(_vec)));

    //std::cout<<(_vec)<<std::endl;
    }

    double shortdis=radius-(radius*cos(GMlib::Angle((M_2PI/bpoints)/2).getRad()));
    //cout<<shortdis<<endl;
    double newradius=radius-(2*shortdis);
    //cout<<newradius<<endl;
    for(int j=1;j<=points;j++)
    {
    auto pt=GMlib::TSVertex<float>(gen_float(-newradius,newradius),gen_float(-newradius,newradius));
    //cout<<pt<<endl;

    if(pt.getParameter().getLength()<newradius)
        {
            bool sign=true;
            for (int k=0; k<this->size(); k++)
                {
                auto op = (*this)[k];
                if ((op.getPos() - pt.getPos()).getLength() < (s/2))//check every points
                                sign = false;
//                if ((this->getElement(k).getParameter() - pt.getParameter().getLength())<s/2)
//                    sign=false;
                }

            if (sign)
                this->insertAlways(pt);
        }
        else
        {
            j--;
        }

    }
      this->triangulateDelaunay();

}

//-------regular points in drum
void regularTriangulate(double beginradius,int cir,int beginpoints){
  //auto p0=GMlib::Point<float,2>(0,0);
    //this->insertAlways(beginradius,0);
  this->insertAlways(GMlib::TSVertex<float>(GMlib::Point<float,2>(0,0)));

  for(int j=1;j<=cir;j++){

    auto _vec=GMlib::Vector<float,2>((float(j)*beginradius),0);


    GMlib::SqMatrix<float,2> nMatrix(GMlib::Angle(M_2PI/(beginpoints*j)),GMlib::Vector<float,2>(1.0f,0.0f),GMlib::Vector<float,2>(0.0f,1.0f));


      for(int i=1;i<=beginpoints*j;i++){
              _vec=(nMatrix *_vec);

              this->insertAlways(GMlib::TSVertex<float>(static_cast<GMlib::Point<float,2>>(_vec)));
      }
      //std::cout<<this->size()<<std::endl;
      //beginpoints+=beginpoints;
}
      this->triangulateDelaunay();
}


//-------transfer type for creat random points

   float gen_float(float a, float b)
   {
       float random = (float)std::rand() / (float) RAND_MAX;
       float diff = b-a;
       float r = random*diff;
       return a+r;
   }

//--------find two triangles have common egde

   GMlib::Vector<GMlib::Vector<float,2>,3> findVector(GMlib::TSEdge<float>* sameEdge){



       auto p0=sameEdge->getFirstVertex();
       auto p1=sameEdge->getLastVertex();
       GMlib::TSVertex<float> *p2;
       GMlib::TSVertex<float> *p3;

       auto triangles=sameEdge->getTriangle();
       auto point2=triangles[0]->getVertices();
       auto point3=triangles[1]->getVertices();

       for(int i=0;i<3;i++){
            if(point2[i]!=p0 &&point2[i]!=p1)
            {
                 p2=point2[i];
            }
            if(point3[i]!=p0 &&point3[i]!=p1)
            {
                 p3=point3[i];
            }
                            }

       GMlib::Vector<float,2> d;
       GMlib::Vector<float,2> a1;
       GMlib::Vector<float,2> a2;

       d=(p1->getParameter())-(p0->getParameter());
       a1=(p2->getParameter())-(p0->getParameter());
       a2=(p3->getParameter())-(p0->getParameter());

       GMlib::Vector<GMlib::Vector<float,2>,3> returnthreeof;
       returnthreeof[0]=d;
       returnthreeof[1]=a1;
       returnthreeof[2]=a2;

       return returnthreeof;
   }

//------find all triangles in same notes

   GMlib::Vector<GMlib::Vector<float,2>,3> findVectorD(GMlib::TSVertex<float> *node,GMlib::TSTriangle<float> *triangle){

           auto tps=triangle->getVertices();

           GMlib::Vector<GMlib::Vector<float,2>,3> returnthreeofD;
           for(int i=0;i<3;i++)
           {

               if(node == tps[1])
               {
                   std::swap(tps[0],tps[1]);
                   std::swap(tps[1],tps[2]);
               }
               if(node == tps[2])
               {
                   std::swap(tps[0],tps[2]);
                   std::swap(tps[1],tps[2]);
               }


           GMlib::Vector<float,2> d1,d2,d3;

           d1=tps[2]->getParameter()-tps[0]->getParameter();
           d2=tps[1]->getParameter()-tps[0]->getParameter();
           d3=tps[2]->getParameter()-tps[1]->getParameter();


           returnthreeofD[0]=d1;
           returnthreeofD[1]=d2;
           returnthreeofD[2]=d3;

           return returnthreeofD;
          }

 }


//-------Calculate the matrix
//-------first initialization matrix, b vector

   void computeM()
   {
       for (int i=0; i<this->size();i++)
       {

           if (!(*this)[i].boundary())
                                                                //if(this->getVertex(i)->boundary())
               _nodes+=&(*this)[i];
                                                                //_nodes += this->getVertex(i);
                                                                //std::cout<<_nodes.size()<<std::endl;

       }
         _Amatrix.setDim(_nodes.size(),_nodes.size());

         for(int i =0;i<_nodes.size();++i){//++i

             for(int j =0;j<_nodes.size();j++)
             {
                 _Amatrix[i][j]=0;
             }
         }

         _b.setDim(_nodes.size());
        for(int i =0;i<_nodes.size();i++){
            _b[i]=0.0f;
        }

       for(int i=0;i<_nodes.size();i++)
        {
           for(int j=i+1; j<_nodes.size();j++)
             {
                   GMlib::TSEdge<float>* sameEdge=nullptr;
                   auto edge1=_nodes[i]->getEdges();
                   auto edge2=_nodes[j]->getEdges();

                  for(int k=0; k<edge1.size();k++)
                  {

                        for(int g=0; g<edge2.size();g++)
                        {
                           if(edge1[k]==edge2[g])
                                {
                                    sameEdge=edge1[k];
                                }
                        }
                   }
              //Computing the stiffness matrix

                   if(sameEdge!=nullptr){
                   auto getVector=findVector(sameEdge);
                   double dd=1/(getVector[0]*getVector[0]);
                   //first trangle
                   double area1=std::abs(getVector[0]^getVector[1]);
                   double dh1=dd*(getVector[1]*getVector[0]);
                   double h1=dd*area1*area1;
                   //sencond triangle
                   double area2=std::abs(getVector[2]^getVector[0]);
                   double dh2=dd*(getVector[2]*getVector[0]);
                   double h2=(dd*area2*area2);
                   //computing Amatrix[i][j]

                   _Amatrix[i][j]=_Amatrix[j][i]=(((dh1*(1-dh1)/h1)-dd)*((area1)/2)+((dh2*(1-dh2)/h2)-dd)*((area2)/2));

                   }
                  else{
                       _Amatrix[i][j]=_Amatrix[j][i]=0;
                   }

              }
           }
       //Computing the stiffness matrix on diagonal

       for(int i=0;i<_nodes.size();i++){
            auto t=_nodes[i]->getTriangles();
                for(int j=0;j<t.size();j++)
                {
                    auto Dvector=findVectorD(_nodes[i],t[j]);
                    _Amatrix[i][i]+=((Dvector[2]*Dvector[2])/(2*std::abs(Dvector[0]^Dvector[1])));
                }
       }

       //Bvector
       for(int i=0; i<_nodes.size();i++){
           auto t=_nodes[i]->getTriangles();
               for(int j=0;j<t.size();j++){
                   auto Bvector=findVectorD(_nodes[i],t[j]);
                   _b[i]+=((std::abs(Bvector[0]^Bvector[1]))/6);
               }
       }


            // ij

             for(int i=0;i<_nodes.size();i++)
             {
                 for(int j=0;j<_nodes.size();j++)
                 {

                    std::cout<<_Amatrix[i][j]<<"  ";
                 }

                 std::cout<< std::endl;
             }

//             cout<<"-------b-------"<<endl;
             for(int i=0;i<_nodes.size();i++){
                 cout<<_b[i]<<endl;
             }

             _MatrixInverted=(_Amatrix.invert());


   }

protected:

//-----creat force for drum

   void localSimulate(double dt)
   {

       if(_swp)
       {
           _force += dt;
           if(_force>1)
               _swp=false;
       }
       else
       {
           _force -= dt;
           if(_force<-1)
               _swp=true;
       }

       _x= _MatrixInverted*_b*_force;
       for(int i=0;i<_nodes.size();i++)
       {
           _nodes[i]->setZ(_x[i]);
       }
       this->replot();

   }


private:


//Array<TSVertex<T>*>     getVertices() const;

GMlib::ArrayLX<GMlib::TSVertex<float>*> _nodes;// 数组
GMlib::DMatrix<float> _Amatrix;
GMlib::DMatrix<float> _MatrixInverted;
GMlib::DVector<float> _b;
//GMlib::TSEdge<float>  _edge;
GMlib::DVector<float> _x;
float _force;
bool _swp;


};


#endif
