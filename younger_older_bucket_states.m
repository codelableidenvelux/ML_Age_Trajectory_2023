%% find buckets

% We onlyconsider the cases when there are 2 stable points and 1 unstable
% point.

pois = (n_stable_points == 2) & (n_unstable_points == 1);
pois_i = find(pois);

% all_mods{pois_i(2)}.plot('potential_eff')

all_buckets = cell(length(pois_i), 2);
pids_for_buckets = sorted_pids(pois_i);
all_dom = cell(length(pois_i), 1);
all_y = cell(length(pois_i), 1);

for ii = 1:length(pois_i)

    IDX = pois_i(ii);

    MU = all_MU(IDX);
    SIGMA = all_SIGMA(IDX);

    mod = all_mods{IDX};
    ueff = mod.potential_eff.ueff;
    U = ueff(mod.potential_eff.dom);
    normU = @(x) ueff(x) - min(U) + 0.1 * (max(U) - min(U));

    mod.equilibria = mod.find_equilibria('effective');
    % extract stable and unstable points
    unstable_eq = mod.equilibria([mod.equilibria.stable] == 0);
    stable_eq = mod.equilibria([mod.equilibria.stable] == 1);

    x_eqL = stable_eq(1).x;
    y_eqL = normU(x_eqL);
    x_eqR = stable_eq(2).x;
    y_eqR = normU(x_eqR);

    x_un = unstable_eq(1).x;
    y_un = normU(x_un);

    yPE = normU(mod.potential_eff.dom);
    yPE = yPE - y_un - 1e-2;
    zcr = find(diff(sign(yPE)) ~= 0);
    delta = mean(abs(diff(mod.potential_eff.dom)));

    [~, x_un_idx] = min(abs(mod.potential_eff.dom - x_un));
    % 3 case
    if length(zcr) == 2
        fprintf('Case 2 zcr\n')
        bucketL = [mod.potential_eff.dom(zcr(1)), x_un - delta] * SIGMA + MU;
        bucketR = [x_un + delta, mod.potential_eff.dom(zcr(2))] * SIGMA + MU;
    elseif length(zcr) == 1
        % here one is the side and one is the zero crossing so we need to check
        % where it is w.r.t. to the x_un
        if (x_un_idx > zcr(1))
            fprintf('Case 1 zcr: RIGHT is side\n')
            yPE = normU(mod.potential_eff.dom);
            yPE = yPE - normU(mod.potential_eff.dom(end-1)) - 1e-10;
            zcr2 = find(diff(sign(yPE)) ~= 0);

            bucketR = [mod.potential_eff.dom(zcr2(end-1)), mod.potential_eff.dom(end-1)] * SIGMA + MU;


            bucketL = [mod.potential_eff.dom(zcr(1)), x_un - delta] * SIGMA + MU;

        else
            fprintf('Case 1 zcr: LEFT is side\n')
            bucketR = [x_un + delta, mod.potential_eff.dom(zcr(1))] * SIGMA + MU;

            yPE = normU(mod.potential_eff.dom);
            yPE = yPE - normU(mod.potential_eff.dom(2)) - 1e-10;
            zcr2 = find(diff(sign(yPE)) ~= 0);
            bucketL = [mod.potential_eff.dom(2), mod.potential_eff.dom(zcr2(2))] * SIGMA + MU;

        end
    else % == 0
        fprintf('Case 0 zcr\n')
        % both zero-crossing are actually the sides
        yPE = normU(mod.potential_eff.dom);
        yPE = yPE - normU(mod.potential_eff.dom(2)) - 1e-10;
        zcr2 = find(diff(sign(yPE)) ~= 0);
        bucketL = [mod.potential_eff.dom(2), mod.potential_eff.dom(zcr2(2))] * SIGMA + MU;

                    yPE = normU(mod.potential_eff.dom);
            yPE = yPE - normU(mod.potential_eff.dom(end-1)) - 1e-10;
            zcr2 = find(diff(sign(yPE)) ~= 0);

            bucketR = [mod.potential_eff.dom(zcr2(end-1)), mod.potential_eff.dom(end-1)] * SIGMA + MU;

    end

    all_buckets{ii, 1} = bucketL;
    all_buckets{ii, 2} = bucketR;
    all_dom{ii, 1} = mod.potential_eff.dom * SIGMA + MU;
    all_y{ii, 1} = normU(mod.potential_eff.dom);

    % subplot(2,1,1)
    %
    % plot(sorted_preds{IDX}, 1:length(sorted_preds{IDX}))
    % hold on
    % plot([bucketL(1), bucketL(1)], [1, length(sorted_preds{IDX})], '--r')
    % plot([bucketL(2), bucketL(2)], [1, length(sorted_preds{IDX})], '--r')
    % plot([bucketR(1), bucketR(1)], [1, length(sorted_preds{IDX})], '--g')
    % plot([bucketR(2), bucketR(2)], [1, length(sorted_preds{IDX})], '--g')
    %
    % subplot(2,1,2)
    % mod.plot('potential_eff')

end
%%
function [slope, width, depth] = calc_slope(x1, y1, x2, y2)
width = x1 - x2;
depth = y1 - y2;
slope = depth / width;
end
