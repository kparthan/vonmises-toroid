function [] = bvm_cosine_norm_constant(kappa1, kappa2, kappa3)

  bess1 = besseli(0,kappa1);
  bess2 = besseli(0,kappa2);
  bess3 = besseli(0,kappa3);
  log_f0 = log(bess1) + log(bess2) + log(bess3);
  sprintf('log_f0: %.6e',log_f0)

  log_even_sum = compute_log_series_sum(2,kappa1,kappa2,kappa3) + log(2);
  sprintf('log_even_sum: %.6e',log_even_sum)

  log_ac = 0;
  if log_f0 > log_even_sum
    ratio = log_even_sum - log_f0;
    tmp = 1 + exp(ratio);
    log_ac = log_f0 + log(tmp);
  else
    ratio = log_f0 - log_even_sum;
    tmp = 1 + exp(ratio);
    log_ac = log_even_sum + log(tmp);
  end
  sprintf('log_ac: %.6e',log_ac)

  %log_odd_sum = compute_log_series_sum(1,kappa1,kappa2,kappa3) + log(2);
  log_odd_sum = compute_log_series_sum_bounded(1,kappa1,kappa2,kappa3,log_ac) + log(2);
  sprintf('log_odd_sum: %.6e',log_odd_sum)

  log_diff = log_odd_sum - log_ac;
  ratio = exp(log_diff);
  sprintf('ratio: %.6e',ratio)  % b/(a+c)

  log_norm_const = 0;
  if kappa3 > 0 % series has all +ve/-ve terms
    tmp = 1 - ratio;
    log_norm_const = (2*log(2*pi)) + log(tmp) + log_ac;
  elseif kappa3 < 0  % series has all +ve terms
    tmp = 1 + ratio;
    log_norm_const = (2*log(2*pi)) + log(tmp) + log_ac;
  end
  sprintf('log_norm_const: %.6e',log_norm_const)

  norm_const = exp(log_norm_const);
  sprintf('norm_const: %.6e',norm_const)

end

function [log_series_sum] = compute_log_series_sum(begin,kappa1,kappa2,kappa3)

  series_sum = 1;
  abs_kappa3 = abs(kappa3);

  j = begin;

  bess1 = besseli(j,kappa1);
  bess2 = besseli(j,kappa2);
  bess3 = besseli(j,abs_kappa3);
  log_f1 = log(bess1) + log(bess2) + log(bess3);

  while 1
    j = j + 2;
    bess1 = besseli(j,kappa1);
    bess2 = besseli(j,kappa2);
    bess3 = besseli(j,abs_kappa3);
    log_fj = log(bess1) + log(bess2) + log(bess3);
    log_tj = log_fj - log_f1;
    tj = exp(log_tj);
    series_sum = series_sum + tj; 
    if (tj / series_sum < eps)
      break;
    end
  end

  log_series_sum = log_f1 + log(series_sum);
  %disp(log_series_sum);

end

function [log_series_sum] = compute_log_series_sum_bounded(begin,kappa1,kappa2,kappa3,log_bound)

  series_sum = 1;
  abs_kappa3 = abs(kappa3);

  j = begin;

  bess1 = besseli(j,kappa1);
  bess2 = besseli(j,kappa2);
  bess3 = besseli(j,abs_kappa3);
  log_f1 = log(bess1) + log(bess2) + log(bess3);

  while 1
    j = j + 2;
    bess1 = besseli(j,kappa1);
    bess2 = besseli(j,kappa2);
    bess3 = besseli(j,abs_kappa3);
    log_fj = log(bess1) + log(bess2) + log(bess3);
    log_tj = log_fj - log_f1;
    tj = exp(log_tj);
    series_sum = series_sum + tj; 
    check = log_f1 + log(series_sum);
    if (check > log_bound)
      series_sum = series_sum - tj;
      break;
    end
    if (tj / series_sum < eps)
      break;
    end
  end

  log_series_sum = log_f1 + log(series_sum);

end
