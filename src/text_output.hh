/*
 * text_output.hh
 *
 *  Created on: 30 Jan 2015
 *      Author: jurak
 */

#ifndef SRC_TEXT_OUTPUT_HH_
#define SRC_TEXT_OUTPUT_HH_

#include <vector>
#include <limits>
#include <string>
#include <stdexcept>

template <typename DGF, typename GV>
class TextOutput{
public:
  TextOutput(DGF const & dgf, GV const & gv) : dgf_(dgf), gv_(gv) { }
  void write(std::string const & filename) const;
private:
  DGF const & dgf_;
  GV  const & gv_;
};

template <typename DGF, typename GV>
void TextOutput<DGF,GV>::write(std::string const & filename) const
{
  const int dim = GV::dimension;
  using Index =  typename GV::IndexSet::IndexType;
  Index max = std::numeric_limits<Index>::max();
  using ctype =  typename GV::ctype;
  Dune::LeafMultipleCodimMultipleGeomTypeMapper<typename GV::Grid, Dune::MCMGVertexLayout> mapper( gv_.grid());
  int nv = mapper.size(); // number of vertices
  std::vector<Index> done(nv,1);
  std::vector<typename DGF::Traits::RangeType>  values(nv);
  std::vector<typename DGF::Traits::DomainType> coordinates(nv);
  for(auto it = gv_.template begin<0>(); it != gv_.template end<0>(); ++it)
  {
      Dune::GeometryType gt = it->type();
      const Dune::template ReferenceElement<ctype,dim> &ref =
                     Dune::ReferenceElements<ctype,dim>::general(gt);
      int vertexsize = it->geometry().corners();
      for(int i=0; i<vertexsize; ++i)
      {
	 typename DGF::Traits::RangeType value = 0;
	 Index vi = mapper.subIndex(*it, i, dim); // vertex global index
	 if(done[vi] == 1){
	     // Index is not yet taken
             auto const& loc_x = ref.position(i, dim); // local coordinate of the vertex
             dgf_.evaluate(*it, loc_x, value);  // function value in the vertex
	     values[vi] = value;
	     coordinates[vi] = it->geometry().corner(i);
	     done[vi]=0; // we are don with this vertex
	 }
	 // do nothing if vertex is already processed
      }
  }
      // Check
      for(auto & x : done) { if(x != 0) throw std::runtime_error("TextOutput<DGF,GV>::write(): internal error.");}

      std::ofstream file(filename);
      for(int i=0; i<nv; ++i){
	  file << coordinates[i] << "   " <<  values[i] << "\n";
//	  std::cout << coordinates[i] << "   " <<  values[i] << "\n";
      }
      file.close();
}

#endif /* SRC_TEXT_OUTPUT_HH_ */
