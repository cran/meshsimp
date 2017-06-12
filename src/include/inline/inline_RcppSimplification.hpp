/*!	\file	inline_RcppSimplification.hpp
	\brief	Implementations of inline methods for class RcppSimplification. */
	
#ifndef HH_INLINERCPPSIMPLIFICATION_HH
#define HH_INLINERCPPSIMPLIFICATION_HH

//
// Get methods
//

INLINE int RcppSimplification::getNumNodes() const
{
	return simplifier.getCPointerToMesh()->getNumNodes();
}


INLINE int RcppSimplification::getNumElems() const
{
	return simplifier.getCPointerToMesh()->getNumElems();
}


INLINE int RcppSimplification::getNumData() const
{
	return simplifier.getCPointerToMesh()->getNumData();
}


//
// Run the simplification
//

INLINE void RcppSimplification::simplify(const UInt & numNodesMax, const CharacterVector & file)
{
	simplifier.simplify(numNodesMax, true, as<string>(file));
}

#endif
