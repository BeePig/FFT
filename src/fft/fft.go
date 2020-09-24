package fft

import (
	//"github.com/mjibson/go-dsp/dsputils"
	//"github.com/mjibson/go-dsp/fft"
	"math"
	"sync"
	"utils"
)

var numcore int = 0
var fft2DChan chan fft1DRow
var transFFT [][]complex128
var w sync.WaitGroup
var rowPerProc int = 0

type fft1DRow struct {
	p     int
	fft1D [][]complex128
}

type data struct {
	x     complex128
	index int
}

var XChan chan data

func FFT(x []float64) (r []complex128) {
	n := len(x)
	r = make([]complex128, n)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			cmpl := complex(x[j], 0)
			r[i] = r[i] + cmpl*getFactor(j*i, n)
		}
	}
	return
}
func BinaryExchange(X []complex128) (Y []complex128) {

	n := len(X)
	Y = make([]complex128, n)
	R := make([]complex128, n)
	r := int(log2(uint(n)))
	mask := make([]int, r)
	for i, _ := range mask {
		mask[i] = getMaskAnd(i, n-1, int(r))
	}
	for i := 0; i < n; i++ {
		R[i] = X[i]
	}
	XChan = make(chan data, n)
	for m := 0; m < int(r); m++ {
		for i := 0; i < n; i++ {
			w.Add(1)
			go binEx(R, i, n, int(r), m, mask[m])
		}
		w.Wait()
		updateFFT(R)
	}
	for i := 0; i < n; i++ {
		revert := revertBit(i, r)
		Y[revert] = R[i]
	}
	return
}

func CooleyTukey(X []complex128) (Y []complex128) {
	n := len(X)
	r := int(log2(uint(n)))
	S := make([]complex128, n)
	R := make([]complex128, n)
	Y = make([]complex128, n)
	for i := 0; i < n; i++ {
		R[i] = X[i]
	}
	for m := 0; m < r; m++ {
		for i := 0; i < n; i++ {
			S[i] = R[i]
		}

		for i := 0; i < n; i++ {
			exp := getExp(i, r, m)
			mask := getMaskAnd(m, n-1, r)
			j := i & mask
			k := j + int(math.Exp2(float64(r-1-m)))
			R[i] = S[j] + S[k]*getFactor(exp, n)
		}
	}

	for i := 0; i < n; i++ {
		revert := revertBit(i, r)
		Y[revert] = R[i]
	}
	return
}

func FFT2DTranspose(x [][]float64) [][]complex128 {
	if numcore == 0 {
		//numcore = runtime.GOMAXPROCS(0)
		numcore = len(x)
	}
	m := len(x)
	fft2DChan = make(chan fft1DRow, numcore)
	if m%numcore != 0 {
		return nil
	}
	transFFT = create2DMatrix(m, m)
	rowPerProc = m / numcore
	for p := 0; p < numcore; p++ {
		w.Add(1)
		go fftRows(p*rowPerProc, rowPerProc, p, utils.ToComplex2(x))
	}
	w.Wait()
	isdone := make(chan bool)
	go fftTranspose(isdone)
	<-isdone
	for p := 0; p < numcore; p++ {
		w.Add(1)
		go fftRows(p*rowPerProc, rowPerProc, p, transFFT)
	}
	w.Wait()
	go fftTranspose(isdone)
	<-isdone
	return transFFT
}

func FFT2D(x [][]float64) (X [][]complex128) {
	n := len(x)

	X = create2DMatrix(n, n)

	fft1D := make([][]complex128, n)
	for row, cols := range x {
		fft1D[row] = BinaryExchange(utils.ToComplex(cols))
	}

	factors := create2DMatrix(n, n)
	for r1 := 0; r1 < n; r1++ {
		for l1 := 0; l1 < n; l1++ {
			factors[r1][l1] = getFactor(r1*l1, n)
		}

	}
	for r1 := 0; r1 < n; r1++ {
		for r2 := 0; r2 < n; r2++ {
			X[r1][r2] = getFFT2DElements(factors, fft1D, n, r1, r2)
		}
	}
	return X
}

func getMaskAnd(index int, maxNum int, maxBit int) int {
	r := maxNum - int(math.Exp2(float64(maxBit-1-index)))
	return r
}

func binEx(X []complex128, i int, n int, r int, m int, mask int) {
	var result complex128 = 0
	var p1 int = i
	j := p1 & mask
	k := j + int(math.Exp2(float64(r-1-m)))

	exp := getExp(i, r, m)
	result = X[j] + X[k]*getFactor(exp, n)
	defer w.Done()
	XChan <- data{result, i}
}

func updateFFT(X []complex128) {
	i := 0
L:
	for {
		select {
		case d := <-XChan:
			X[d.index] = d.x
			i++
		default:
			break L

		}
	}
}

func getExp(i int, r int, m int) int {
	if m == r-1 {
		return revertBit(i, r)
	}
	maxNum := int(math.Exp2(float64(r))) - 1
	mask := getMaskAnd(m+1, maxNum, r)
	revert := revertBit(i&mask, r)
	return int(revert << (r - 1 - m))
}

func fftRows(startRow int, q int, p int, x [][]complex128) {
	defer w.Done()
	ffts := create2DMatrix(q, len(x))
	i := 0
	for r := startRow; r < startRow+q; r++ {
		ffts[i] = CooleyTukey(x[r])
		i++
	}
	fft2DChan <- fft1DRow{
		p:     p,
		fft1D: ffts,
	}
}

func fftTranspose(isdone chan bool) {
L:
	for {
		select {
		case fft := <-fft2DChan:
			contributeFFTTranspose(fft, transFFT)
		default:
			break L
		}
	}
	l := len(transFFT)
	temp := create2DMatrix(l, l)
	for i := 0; i < l; i++ {
		for j := 0; j < l; j++ {
			temp[i][j] = transFFT[j][i]
		}
	}
	copy(transFFT, temp)
	isdone <- true
}

func transDistributedFFT() {
L:
	for {
		select {
		case fft := <-fft2DChan:
			contributeFFTTranspose(fft, transFFT)
		default:
			break L
		}
	}
}
func contributeFFTTranspose(data fft1DRow, transFFT [][]complex128) {
	index := 0
	for i := 0; i < data.p; i++ {
		index = index + rowPerProc
	}
	for i := 0; i < rowPerProc; i++ {
		transFFT[index] = data.fft1D[i]
		index++
	}
}

func getFFT2DElements(factors [][]complex128, fft1D [][]complex128, cols int, r1 int, r2 int) complex128 {
	var x complex128
	for i := 0; i < cols; i++ {
		x = x + factors[r1][i]*fft1D[i][r2]
	}
	return x
}
func create2DMatrix(m int, n int) [][]complex128 {
	r := make([][]complex128, m)
	for i, _ := range r {
		r[i] = make([]complex128, n)
	}
	return r
}

func getFactor(exp int, n int) complex128 {
	sin, cos := math.Sincos(-2 * math.Pi / float64(n) * float64(exp))
	return complex(math.Round(cos*100000)/100000, math.Round(sin*100000)/100000)
}

// log2 returns the log base 2 of v
// from: http://graphics.stanford.edu/~seander/bithacks.html#IntegerLogObvious
func log2(v uint) uint {
	var r uint

	for v >>= 1; v != 0; v >>= 1 {
		r++
	}

	return r
}

func revertBit(n int, maxBit int) (result int) {
	masks := make([]int, maxBit)
	a := 1
	temp := make([]int, maxBit)
	for i := 0; i < maxBit; i++ {
		for j := 0; j < maxBit-1-i; j++ {
			a = a << 1
		}
		masks[i] = a
		a = 1
	}

	for i := 0; i < maxBit; i++ {
		temp[i] = (masks[i] & n) >> (maxBit - 1 - i)
	}

	for i, _ := range temp {
		result = result + temp[i]*int(math.Exp2(float64(i)))
	}
	return
}
