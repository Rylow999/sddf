// =====================================================================
// Calculadora de Curvatura Espectral SDDF  -- version corregida
//
//     G[u] = INT [ d(ln E)/d(ln k) ]^2  d(ln k)
//
// Cambios respecto de la version original:
//   1) Se elimina el piso absoluto `1e-30` del lector, que aplastaba la
//      cola de disipacion en una meseta plana y generaba un salto espurio
//      en la derivada logaritmica.
//   2) Se agregan dos truncados de dominio, porque la cola de disipacion
//      del espectro de Pope tiene pendiente log-log DIVERGENTE
//      (s(k) = -5/3 - (4*beta/3)*(k*eta)^(4/3)), asi que la integral no
//      converge al extender k: hay que definir explicitamente la ventana.
//        - `truncate_tail`      : corte por amplitud relativa (rapido)
//        - `truncate_by_slope`  : corte por desviacion de pendiente
//                                 respecto de -5/3  [RECOMENDADO]
//   3) Se acepta el archivo por argumento de linea de comandos.
//
// Uso:  cargo run --release -- spectrum_Re_100.csv [delta]
// =====================================================================

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

const SLOPE_K41: f64 = -5.0 / 3.0;

/// Espectro sintetico de Pao (equivalente funcional a Pope con f_eta
/// exponencial simple). beta=2.25 (=1.5*C_K con C_K=1.5): fallback
/// quimera acordado en v3.1. NO usar 5.2 (Pope) sin verificar el espectro.
fn pope_spectrum(k: f64, epsilon: f64, nu: f64, c: f64, beta: f64) -> f64 {
    let eta: f64 = (nu.powi(3) / epsilon).powf(0.25);
    let k_eta: f64 = k * eta;
    let f_eta: f64 = (-beta * k_eta.powf(4.0 / 3.0)).exp();
    let e_inertial: f64 = c * epsilon.powf(2.0 / 3.0) * k.powf(-5.0 / 3.0);
    e_inertial * f_eta
}

/// Pendiente log-log local: diferencias centradas en el interior,
/// forward/backward en los bordes. Identica a la version original.
fn log_slope(k_vec: &[f64], e_vec: &[f64]) -> Vec<f64> {
    let n: usize = k_vec.len();
    let mut deriv: Vec<f64> = vec![0.0; n];
    if n < 2 {
        return deriv;
    }

    let log_k: Vec<f64> = k_vec.iter().map(|x| x.ln()).collect();
    let log_e: Vec<f64> = e_vec.iter().map(|x| x.ln()).collect();

    for i in 1..(n - 1) {
        let dx: f64 = log_k[i + 1] - log_k[i - 1];
        let dy: f64 = log_e[i + 1] - log_e[i - 1];
        deriv[i] = if dx.abs() > 1e-12 { dy / dx } else { 0.0 };
    }

    deriv[0] = (log_e[1] - log_e[0]) / (log_k[1] - log_k[0]);
    deriv[n - 1] = (log_e[n - 1] - log_e[n - 2]) / (log_k[n - 1] - log_k[n - 2]);

    deriv
}

/// G[u] por regla del trapecio sobre ln k. Identica a la original.
fn compute_spectral_curvature(k_vec: &[f64], e_vec: &[f64]) -> f64 {
    let n: usize = k_vec.len();
    if n < 3 {
        return 0.0;
    }

    let log_k: Vec<f64> = k_vec.iter().map(|x| x.ln()).collect();
    let deriv: Vec<f64> = log_slope(k_vec, e_vec);

    let mut g_total: f64 = 0.0;
    for i in 0..(n - 1) {
        let dx: f64 = log_k[i + 1] - log_k[i];
        let y1: f64 = deriv[i].powi(2);
        let y2: f64 = deriv[i + 1].powi(2);
        g_total += 0.5 * (y1 + y2) * dx;
    }

    g_total
}

/// Corte por amplitud relativa: recorta en el primer k donde
/// E(k) < rel_thresh * E(k_0). No necesita conocer nu ni epsilon.
/// OJO: el exponente que se ajuste despues DEPENDE de rel_thresh.
fn truncate_tail(k_vec: &[f64], e_vec: &[f64], rel_thresh: f64) -> (Vec<f64>, Vec<f64>) {
    let n: usize = k_vec.len();
    if n < 3 {
        return (k_vec.to_vec(), e_vec.to_vec());
    }
    let thresh: f64 = rel_thresh * e_vec[0];
    let mut cut: usize = n;
    for i in 0..n {
        if e_vec[i] < thresh {
            cut = i;
            break;
        }
    }
    if cut < 3 {
        cut = 3;
    }
    (k_vec[..cut].to_vec(), e_vec[..cut].to_vec())
}

/// Corte por desviacion de pendiente [RECOMENDADO]: recorta en el primer k
/// donde |d(ln E)/d(ln k) + 5/3| > delta. Es data-driven y equivale a un
/// corte en k*eta = cte, que es lo fisicamente correcto: aisla el rango
/// inercial y deja fuera la zona de disipacion, cuya contribucion diverge.
/// Corte por desviacion de pendiente CON INTERPOLACION LINEAL en (ln k, s),
/// sincronizado con sddf_core.truncate_by_slope_interp (v3.1).
/// El corte se hace contra `s_ref` (K41 por defecto), no contra el q que
/// se este estimando. Devuelve (k_trunc, e_trunc, n_eliminados, ok, k_corte).
fn truncate_by_slope_interp(
    k_vec: &[f64],
    e_vec: &[f64],
    delta: f64,
    s_ref: f64,
) -> (Vec<f64>, Vec<f64>, usize, bool, f64) {
    let n: usize = k_vec.len();
    if n < 3 {
        return (k_vec.to_vec(), e_vec.to_vec(), 0, false, f64::NAN);
    }
    let s: Vec<f64> = log_slope(k_vec, e_vec);
    let dev: Vec<f64> = s.iter().map(|x| (x - s_ref).abs()).collect();
    let j_opt: Option<usize> = dev.iter().position(|&d| d > delta);
    let j: usize = match j_opt {
        None => return (k_vec.to_vec(), e_vec.to_vec(), 0, true, *k_vec.last().unwrap()),
        Some(j) => j,
    };
    if j < 3 {
        return (k_vec[..3].to_vec(), e_vec[..3].to_vec(), n - 3, false, k_vec[2]);
    }
    let lk: Vec<f64> = k_vec.iter().map(|x| x.ln()).collect();
    let d0: f64 = dev[j - 1];
    let d1: f64 = dev[j];
    let mut w: f64 = if d1 == d0 { 0.0 } else { (delta - d0) / (d1 - d0) };
    w = w.clamp(0.0, 1.0);
    let lkc: f64 = lk[j - 1] + w * (lk[j] - lk[j - 1]);
    let lec: f64 = e_vec[j - 1].ln() + w * (e_vec[j].ln() - e_vec[j - 1].ln());
    let mut k_new: Vec<f64> = k_vec[..j].to_vec();
    let mut e_new: Vec<f64> = e_vec[..j].to_vec();
    k_new.push(lkc.exp());
    e_new.push(lec.exp());
    let kc: f64 = lkc.exp();
    (k_new, e_new, n - j, true, kc)
}

/// Lector sin piso absoluto: descarta puntos no finitos o no positivos
/// (el logaritmo no esta definido ahi) en lugar de aplastarlos.
fn read_spectrum_csv(filename: &str) -> Option<(Vec<f64>, Vec<f64>, usize)> {
    let file = File::open(filename).ok()?;
    let reader = BufReader::new(file);
    let mut k_vec: Vec<f64> = Vec::new();
    let mut e_vec: Vec<f64> = Vec::new();
    let mut dropped: usize = 0;

    for line in reader.lines() {
        let line = line.ok()?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = trimmed.split(',').collect();
        if parts.len() >= 2 {
            let k_res: Result<f64, _> = parts[0].trim().parse();
            let e_res: Result<f64, _> = parts[1].trim().parse();
            if let (Ok(k), Ok(e)) = (k_res, e_res) {
                // Sin piso: E<=0 o no finito no se aplasta, se descarta.
                if k.is_finite() && k > 0.0 && e.is_finite() && e > 0.0 {
                    k_vec.push(k);
                    e_vec.push(e);
                } else {
                    dropped += 1;
                }
            }
        }
    }

    if k_vec.len() > 2 {
        Some((k_vec, e_vec, dropped))
    } else {
        None
    }
}

fn main() {
    println!("=== Calculadora de Curvatura Espectral SDDF (corregida) ===");

    let args: Vec<String> = env::args().collect();
    let csv_file: String = if args.len() > 1 {
        args[1].clone()
    } else {
        String::from("spectrum.csv")
    };
    let delta: f64 = if args.len() > 2 {
        args[2].parse().unwrap_or(0.5)
    } else {
        0.5
    };
    let rel_thresh: f64 = 1e-6;

    let (k_vec, e_vec, source) = match read_spectrum_csv(&csv_file) {
        Some((k, e, dropped)) => {
            println!("[OK] Datos cargados desde '{}'. Puntos: {}", csv_file, k.len());
            if dropped > 0 {
                println!("[WARN] {} puntos descartados (E<=0 o no finito).", dropped);
            }
            (k, e, "Real (CSV)")
        }
        None => {
            println!("[INFO] No se encontro '{}'. Generando espectro sintetico de Pope...", csv_file);
            let n_points: usize = 2000;
            let mut k: Vec<f64> = vec![0.0; n_points];
            let mut e: Vec<f64> = vec![0.0; n_points];

            let k_min: f64 = 1.0;
            let k_max: f64 = 1000.0;

            for i in 0..n_points {
                let ratio: f64 = (k_max / k_min).powf(i as f64 / (n_points as f64 - 1.0));
                k[i] = k_min * ratio;
                e[i] = pope_spectrum(k[i], 1.0, 0.01, 1.5, 2.25);
            }
            (k, e, "Sintetico (Pope)")
        }
    };

    // --- Sin truncar: diverge si hay cola de disipacion resuelta --------
    let g_raw: f64 = compute_spectral_curvature(&k_vec, &e_vec);

    // --- Corte por amplitud relativa ------------------------------------
    let (k_rel, e_rel) = truncate_tail(&k_vec, &e_vec, rel_thresh);
    let g_rel: f64 = compute_spectral_curvature(&k_rel, &e_rel);

    // --- Corte por pendiente con interpolacion [recomendado] ------------
    let (k_sl, e_sl, _n_elim, ok_sl, k_corte) =
        truncate_by_slope_interp(&k_vec, &e_vec, delta, SLOPE_K41);
    let g_sl: f64 = compute_spectral_curvature(&k_sl, &e_sl);

    println!("Fuente               : {}", source);
    println!(
        "Rango k completo     : {:.3e} a {:.3e}  ({} pts)",
        k_vec.first().unwrap(),
        k_vec.last().unwrap(),
        k_vec.len()
    );
    println!("G[u] sin truncar     : {:.4}   <-- depende de k_max, no usar", g_raw);
    println!(
        "G*[u] corte amplitud : {:.4}   (thresh={:.0e}, k_max={:.3e}, {} pts)",
        g_rel,
        rel_thresh,
        k_rel.last().unwrap(),
        k_rel.len()
    );
    println!(
        "G*[u] corte pendiente: {:.4}   (delta={:.2}, k_max={:.3e}, {} pts, interp, ok={})  [RECOMENDADO]",
        g_sl,
        delta,
        k_corte,
        k_sl.len(),
        ok_sl
    );

    if let Ok(mut file) = File::create("sddf_history.log") {
        let _ = writeln!(
            file,
            "Source: {}, File: {}, Points: {}, G_raw: {:.4}, G_rel: {:.4}, G_slope: {:.4}, delta: {:.2}",
            source, csv_file, k_vec.len(), g_raw, g_rel, g_sl, delta
        );
    }
}
