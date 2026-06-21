use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Debug)]
pub struct OpacityTable {
    pub xs: Vec<f64>,
    pub log_ts: Vec<f64>,
    pub log_rs: Vec<f64>,
    // 3D grid: [x_idx][log_t_idx][log_r_idx]
    pub kappas: Vec<Vec<Vec<f64>>>,
}

impl OpacityTable {
    pub fn load_from_file<P: AsRef<Path>>(path: P) -> std::io::Result<Self> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);

        let mut xs = Vec::new();
        let mut log_ts = Vec::new();
        let mut log_rs = Vec::new();
        let mut kappas = Vec::new();

        let mut current_x_idx = 0;
        let mut in_table = false;
        let mut parsed_log_rs = false;
        
        let mut current_kappas = Vec::new();

        for line in reader.lines() {
            let line = line?;
            let line = line.trim();

            if line.starts_with("Table for X =") {
                if in_table && !current_kappas.is_empty() {
                    kappas.push(current_kappas.clone());
                    current_kappas.clear();
                    current_x_idx += 1;
                }
                
                let parts: Vec<&str> = line.split('=').collect();
                if parts.len() == 2 {
                    let x_val = parts[1].trim().parse::<f64>().unwrap();
                    xs.push(x_val);
                }
                in_table = true;
                continue;
            }

            if in_table {
                if line.starts_with("logT / logR |") {
                    if !parsed_log_rs {
                        let parts: Vec<&str> = line.split('|').collect();
                        if parts.len() == 2 {
                            for val in parts[1].split_whitespace() {
                                log_rs.push(val.parse::<f64>().unwrap());
                            }
                            parsed_log_rs = true;
                        }
                    }
                    continue;
                }

                if line.starts_with("---") {
                    continue;
                }

                if line.is_empty() {
                    if !current_kappas.is_empty() {
                        kappas.push(current_kappas.clone());
                        current_kappas.clear();
                        current_x_idx += 1;
                        in_table = false;
                    }
                    continue;
                }

                if line.contains('|') {
                    let parts: Vec<&str> = line.split('|').collect();
                    if parts.len() == 2 {
                        let log_t = parts[0].trim().parse::<f64>().unwrap();
                        if current_x_idx == 0 {
                            log_ts.push(log_t);
                        }
                        
                        let mut row = Vec::new();
                        for val in parts[1].split_whitespace() {
                            row.push(val.parse::<f64>().unwrap());
                        }
                        current_kappas.push(row);
                    }
                }
            }
        }
        
        // Push the last table if the file doesn't end with a blank line
        if in_table && !current_kappas.is_empty() {
            kappas.push(current_kappas);
        }

        Ok(OpacityTable {
            xs,
            log_ts,
            log_rs,
            kappas,
        })
    }

    // Trilinear interpolation for log(kappa) given X, log(T), and log(R)
    pub fn get_log_kappa(&self, x: f64, log_t: f64, log_r: f64) -> f64 {
        // Clamp inputs to table bounds
        let x = x.clamp(*self.xs.first().unwrap(), *self.xs.last().unwrap());
        let log_t = log_t.clamp(*self.log_ts.first().unwrap(), *self.log_ts.last().unwrap());
        let log_r = log_r.clamp(*self.log_rs.first().unwrap(), *self.log_rs.last().unwrap());

        let (ix0, ix1, tx) = Self::find_indices_and_weight(&self.xs, x);
        let (it0, it1, tt) = Self::find_indices_and_weight(&self.log_ts, log_t);
        let (ir0, ir1, tr) = Self::find_indices_and_weight(&self.log_rs, log_r);

        let v000 = self.kappas[ix0][it0][ir0];
        let v001 = self.kappas[ix0][it0][ir1];
        let v010 = self.kappas[ix0][it1][ir0];
        let v011 = self.kappas[ix0][it1][ir1];
        let v100 = self.kappas[ix1][it0][ir0];
        let v101 = self.kappas[ix1][it0][ir1];
        let v110 = self.kappas[ix1][it1][ir0];
        let v111 = self.kappas[ix1][it1][ir1];

        // Interpolate over R
        let c00 = v000 * (1.0 - tr) + v001 * tr;
        let c01 = v010 * (1.0 - tr) + v011 * tr;
        let c10 = v100 * (1.0 - tr) + v101 * tr;
        let c11 = v110 * (1.0 - tr) + v111 * tr;

        // Interpolate over T
        let c0 = c00 * (1.0 - tt) + c01 * tt;
        let c1 = c10 * (1.0 - tt) + c11 * tt;

        // Interpolate over X
        c0 * (1.0 - tx) + c1 * tx
    }

    fn find_indices_and_weight(arr: &[f64], val: f64) -> (usize, usize, f64) {
        if val <= arr[0] {
            return (0, 0, 0.0);
        }
        let last_idx = arr.len() - 1;
        if val >= arr[last_idx] {
            return (last_idx, last_idx, 0.0);
        }

        match arr.binary_search_by(|a| a.partial_cmp(&val).unwrap()) {
            Ok(idx) => (idx, idx, 0.0),
            Err(idx) => {
                let i0 = idx - 1;
                let i1 = idx;
                let t = (val - arr[i0]) / (arr[i1] - arr[i0]);
                (i0, i1, t)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_indices() {
        let arr = vec![1.0, 2.0, 3.0, 4.0];
        let (i0, i1, t) = OpacityTable::find_indices_and_weight(&arr, 2.5);
        assert_eq!(i0, 1);
        assert_eq!(i1, 2);
        assert!((t - 0.5).abs() < 1e-6);
    }
}
