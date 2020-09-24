package utils


func ToComplex(x []float64) []complex128 {
	y := make([]complex128, len(x))
	for n, v := range x {
		y[n] = complex(v, 0)
	}
	return y
}

func ToComplex2(x [][]float64) [][]complex128 {
	y := make([][]complex128, len(x))
	for n, v := range x {
		y[n] = ToComplex(v)
	}
	return y
}
