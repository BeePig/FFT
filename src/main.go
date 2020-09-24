package main

import (
	"fft"
	"fmt"
	"utils"
	"math/rand"
	"time"
)

const LEN_FFT int = 30000000
const LEN_FFT_ROW int = 5477
const MAX_NUM float64 = 10.0

func main() {
	compareCooleyTukey_BinaryExchange()
	compareFFT2D_2DTranspose()
}

func compareCooleyTukey_BinaryExchange(){
	a := make([]float64,LEN_FFT)
	for i:=0;i<LEN_FFT;i++{
		a[i] = rand.Float64()*MAX_NUM
	}
	timeStart := time.Now()
	fft.CooleyTukey(utils.ToComplex(a))
	elapse := time.Since(timeStart).Seconds()
	fmt.Println("Cooley-Tukey run in: ", elapse,"s")
	timeStart = time.Now()
	fft.BinaryExchange(utils.ToComplex(a))
	elapse = time.Since(timeStart).Seconds()
	fmt.Println("BinaryExchange run in: ", elapse,"s")
}

func compareFFT2D_2DTranspose(){
	a := make([][]float64,LEN_FFT_ROW)
	for i:=0;i<LEN_FFT_ROW;i++{
		a[i] = make([]float64, LEN_FFT_ROW)
		for j := 0;j<LEN_FFT_ROW;j++{
			a[i][j] = rand.Float64()*MAX_NUM
		}

	}
	timeStart := time.Now()
	fft.FFT2D(a)
	elapse := time.Since(timeStart).Seconds()
	fmt.Println("FFT2D run in: ", elapse,"s")
	timeStart = time.Now()
	fft.FFT2DTranspose(a)
	elapse = time.Since(timeStart).Seconds()
	fmt.Println("FFT2DTranspose run in: ", elapse,"s")
}


