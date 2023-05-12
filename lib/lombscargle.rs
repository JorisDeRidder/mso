
use std::f64::consts::PI as PI_f64;


/// Compute a Lomb-Scargle frequencygramme.
///
/// myDescription 
/// $$ E = m c^2 $$
/// 
/// # Arguments
///
/// * `time`          - An array with time points, not necessarily equidistant. 
///                     No reference time t0 is subtracted; you need to do this yourself before you call this function.
/// * `signal`        - An array with the corresponding signal points. The mean does not need to be already subtracted. 
/// * `weights`       - An array with the weights that are normalized to sum to 1. $w_n = 1/(W * sigma_n^2)$ with $W = \sum 1/sigma_n^2$. 
/// * `freq_start`    - The first frequency > 0, for which the spectrum needs to be evaluated. Not angular frequency.
///                     Units: reciprocal unit of `time`, i.e. cycles per unit of `time`. Must be strictly larger than 0.0.
/// * `freq_step`     - The step size of the equidistant frequency grid of the frequencygramme. Units: same as `freq_start`. 
/// * `num_freq`      - The number of points in the frequency grid.
/// * `with_constant` - If `true`, model the signal as sine wave + constant, if `false`, only use a sine wave. 
///
/// # Output
///
/// * `spectrum`      - Vector of length `num_freq` containing $(\chi^2_0 - \chi^2(\omega))/\chi^2_0$.
/// * `amplitude_cos` - Vector of length `num_freq` containing the fitted amplitude of the cosine term.
/// * `amplitude_sin` - Vector of length `num_freq` containing the fitted amplitude of the sine term.
/// * `constant`      - Vector of length `num_freq` containing the fitted value for the with_constant
///                     If `with_constant == false` then this is a constant array containing the weighted mean of the y.
///
/// # Panics
/// 
/// -  Panics when time or signal have zero length, when the frequency step is <= 0, or when num_freq <= 0.
///
/// # References
///
/// - Zechmeister & Kürster, 2009, A&A 496, p. 577Z           (This function implements Eqs (4)-(15))
/// - Reegen, 2007, A&A 467, p. 1353
/// - Lomb, 1976, Ap&SS 39, p. 447
/// - Scargle, 1982, ApJ 263, p. 835
///
///
pub fn lombscargle(time: &[f64], signal: &[f64], weights: &[f64], freq_start: f64, freq_step: f64, num_freq: usize, with_constant: bool) 
    -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {

    assert!(time.len() > 0);
    assert!(signal.len() == time.len());
    assert!(freq_step > 0.0);
    assert!(num_freq >= 1);

    let mean = signal.iter().zip(weights.iter()).map(|(x,y)| x*y).sum::<f64>(); // weights must be normalized
    let y: Vec<f64> = signal.iter().map(|x| x-mean).collect();                  // Subtract the weighted mean

    let omega_start = freq_start * 2.0 * PI_f64;                                // ω = 2πν
    let omega_step  = freq_step  * 2.0 * PI_f64;                                // δω = 2π δν

    let mut sum_y        = 0.0;                                                 // Σ w y
    let mut sum_yy       = 0.0;                                                 // Σ w y^2
    let mut sum_ysinx    = vec![0.0; num_freq];                                 // Σ w sin(ωt) y
    let mut sum_ycosx    = vec![0.0; num_freq];                                 // Σ w cos(ωt) y
    let mut sum_sinx     = vec![0.0; num_freq];                                 // Σ w sin(ωt)
    let mut sum_cosx     = vec![0.0; num_freq];                                 // Σ w cos(ωt)
    let mut sum_sinxsinx = vec![0.0; num_freq];                                 // Σ w sin(ωt)^2
    let mut sum_sinxcosx = vec![0.0; num_freq];                                 // Σ w sin(ωt) cos(ωt)

    for n in 0..time.len() {
        sum_y  += weights[n] * y[n];                                            // Should be 0 if weighted mean was subtracted from signal
        sum_yy += weights[n] * y[n] * y[n];
        let sindx = (omega_step*time[n]).sin();                                 // sin(δω t)
        let cosdx = (omega_step*time[n]).cos();                                 // cos(δω t)
        let mut sinx = (omega_start*time[n]).sin();                             // sin(ωt)
        let mut cosx = (omega_start*time[n]).cos();                             // cos(ωt)

        for j in 0..num_freq {
            sum_ysinx[j]    += weights[n] * y[n] * sinx;
            sum_ycosx[j]    += weights[n] * y[n] * cosx;
            sum_sinx[j]     += weights[n] * sinx;            
            sum_cosx[j]     += weights[n] * cosx;            
            sum_sinxsinx[j] += weights[n] * sinx * sinx;     
            sum_sinxcosx[j] += weights[n] * sinx * cosx;     

            // Update the values for sinx and cosx for the _next_ iteration.
            // Use update equations that don't require sin() or cos() functions, 
            // except for every now and then to counter numerical inaccuracies.

            if j % 5000  == 0 {                                                       
                sinx = (time[n]*(omega_start+(j as f64 + 1.0)*omega_step)).sin();
                cosx = (time[n]*(omega_start+(j as f64 + 1.0)*omega_step)).cos();

            } else {
                let tmp = cosx;
                cosx = tmp * cosdx - sinx * sindx;
                sinx = sinx * cosdx + tmp * sindx;
            }
        }
    }

    let mut spectrum:      Vec<f64> = Vec::with_capacity(num_freq);
    let mut amplitude_sin: Vec<f64> = Vec::with_capacity(num_freq);
    let mut amplitude_cos: Vec<f64> = Vec::with_capacity(num_freq);
    let mut constant:      Vec<f64> = Vec::with_capacity(num_freq);

    let yy = sum_yy - sum_y * sum_y;
    for j in 0..num_freq {
        let ys = sum_ysinx[j] - sum_y * sum_sinx[j];
        let yc = sum_ycosx[j] - sum_y * sum_cosx[j];
        if with_constant {
            let ss = sum_sinxsinx[j] - sum_sinx[j] * sum_sinx[j];
            let cc = (1.0 - sum_sinxsinx[j]) - sum_cosx[j] * sum_cosx[j];
            let cs = sum_sinxcosx[j] - sum_cosx[j] * sum_sinx[j];
            let d = cc*ss-cs*cs;
            spectrum.push( (ss*yc*yc + cc*ys*ys - 2.0*cs*yc*ys) / (yy*d) );
            let ampl_cos = (yc*ss-ys*cs) / d;                                                // Eq. (A4) in Zechmeister & Kuerster 
            let ampl_sin = (ys*cc-yc*cs) / d;                                                // Idem.
            amplitude_cos.push(ampl_cos);
            amplitude_sin.push(ampl_sin);
            constant.push(sum_y - ampl_cos * sum_cosx[j] - ampl_sin * sum_sinx[j] + mean);
        } else {
            let ss = sum_sinxsinx[j];
            let cc = 1.0 - sum_sinxsinx[j];
            let cs = sum_sinxcosx[j];
            let d = cc*ss-cs*cs;
            spectrum.push( (ss*yc*yc + cc*ys*ys - 2.0*cs*yc*ys) / (yy*d) );
            amplitude_cos.push((yc*ss-ys*cs) / d);
            amplitude_sin.push((ys*cc-yc*cs) / d);
            constant.push(mean);                                                    // A constant array, for consistency
        }
    }

    (spectrum, amplitude_cos, amplitude_sin, constant)
}



