#include <type_traits>
#include <vector>

namespace PIC_PLUS_PLUS {
template <typename VecType> class Flat2DHelper {
private:
  std::vector<VecType> &m_vec;
  size_t m_curRow;
  size_t m_numRows;
  size_t m_numColumns;

public:
  Flat2DHelper(std::vector<VecType> &vec, size_t numRows, size_t numColumns,
               size_t curRow)
      : m_vec(vec), m_numRows(numRows), m_numColumns(numColumns),
        m_curRow(curRow) {}
  template <typename IntType>
  std::enable_if_t<std::is_integral_v<IntType>, VecType> &
  operator[](IntType column) {
    return m_vec[m_curRow * m_numColumns + column];
  }
  ~Flat2DHelper() = default;
};

template <typename T> class Flat2DVector {
private:
  size_t m_numRows;
  size_t m_numColumns;

public:
  Flat2DVector(size_t rows, size_t columns)
      : m_numRows(rows), m_numColumns(columns) {
    flatVector.reserve(rows * columns);
  }
  void push_back(T &element);
  template <class... Args> void emplace_back(Args &&...args) {
    flatVector.emplace_back(std::forward<Args>(args)...);
  }

  template <typename I>
  std::enable_if_t<std::is_integral_v<I>> reserve(I amountToReserve) {
    flatVector.reserve(amountToReserve);
  }

  template <typename IntType>
  std::enable_if_t<std::is_integral_v<IntType>, Flat2DHelper<T>>
  operator[](IntType row) {
    return Flat2DHelper<T>(flatVector, m_numRows, m_numColumns, row);
  }

  ~Flat2DVector() = default;

private:
  std::vector<T> flatVector;
};
} // namespace PIC_PLUS_PLUS
