#pragma once

template <typename _RawArray, size_t _Dim>
class ArrayView {
public:
typedef _RawArray RawArray;
static const size_t Dim = _Dim;

	explicit ArrayView(RawArray& array...): _data(array) {
		va_list args;
		va_start(args, array);
		
		size_t size = 1;
		for (size_t d = 0 ; d < Dim ; ++d) {
			size_t n = va_arg(args, size_t);
			_shape[d] = n;
			size *= n;
		}
		
		_bound[Dim] = 1;
		for (size_t d = Dim-1 ; d < Dim ; --d) {
			_bound[d] = _bound[d+1] * _shape[d];
		}
	}
	
	Scalar operator()(
	
private:
	unsigned int _flag;
	RawArray& _data;
	std::array<size_t, Dim> _shape;
	std::array<size_t, Dim+1> _bound;
};
