
namespace debug {

template <size_t Rank>
void
debugshow(const char* name,
          const array::ArrayStructure<Rank>& arr,
          std::ostream& os,
          size_t indent_level = 0)
{
  for (size_t i = 0 ; i < indent_level ; ++i) {
    os << "  ";
  }
  os << name << "[" << typeid(arr).name() << "]" << std::endl;
  debugshow("length", arr.length(), os, indent_level+1);
  debugshow("shape",  arr.shape(), os, indent_level+1);
  debugshow("stride", arr.stride(), os, indent_level+1);
}

template <typename ValueType, size_t Rank>
void
debugshow(const char* name,
          const array::Array<ValueType, Rank>& arr,
          std::ostream& os,
          size_t indent_level = 0)
{
  for (size_t i = 0 ; i < indent_level ; ++i) {
    os << "  ";
  }
  os << name << "[" << typeid(arr).name() << "]" << std::endl;
  debugshow("structure", arr.structure(), os, indent_level+1);
  debugshow("data",  arr.data(), os, indent_level+1);

  for (size_t i = 0 ; i < indent_level+1 ; ++i) { os << "  "; }
  os << "*data" << std::endl;

  for (size_t i = 0 ; i < arr.length() ; ++i) {
    for (size_t k = 0 ; k < indent_level+2 ; ++k) {
      os << "  ";
    }
    os << arr.data().get()[i] << std::endl;
  }
  debugshow("own", arr.own(), os, indent_level+1);
}


} // namespace debug
