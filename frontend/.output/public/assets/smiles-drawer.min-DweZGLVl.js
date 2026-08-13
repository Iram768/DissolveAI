var e = Object.create,
  t = Object.defineProperty,
  n = Object.getOwnPropertyDescriptor,
  r = Object.getOwnPropertyNames,
  i = Object.getPrototypeOf,
  a = Object.prototype.hasOwnProperty,
  o = (e, n, r) =>
    n in e
      ? t(e, n, { enumerable: !0, configurable: !0, writable: !0, value: r })
      : (e[n] = r),
  s = (e, t) => () => {
    try {
      return (t || e((t = { exports: {} }).exports, t), t.exports);
    } catch (e) {
      throw ((t = 0), e);
    }
  },
  c = (e, i, o, s) => {
    if ((i && typeof i == `object`) || typeof i == `function`)
      for (let c of r(i))
        !a.call(e, c) &&
          c !== o &&
          t(e, c, {
            get: () => i[c],
            enumerable: !(s = n(i, c)) || s.enumerable,
          });
    return e;
  },
  l = (n, r, a) => (
    (a = n == null ? {} : e(i(n))),
    c(
      r || !n || !n.__esModule
        ? t(a, `default`, { value: n, enumerable: !0 })
        : a,
      n,
    )
  ),
  u = (e, t, n) => o(e, typeof t == `symbol` ? t : t + ``, n),
  d = s((e, t) => {
    (function (n, r) {
      typeof e == `object` && typeof t < `u`
        ? (t.exports = r())
        : typeof define == `function` && define.amd
          ? define(r)
          : ((n = typeof globalThis < `u` ? globalThis : n || self),
            (n.chroma = r()));
    })(e, function () {
      for (
        var e = function (e, t, n) {
            return (
              t === void 0 && (t = 0),
              n === void 0 && (n = 1),
              e < t ? t : e > n ? n : e
            );
          },
          t = e,
          n = function (e) {
            ((e._clipped = !1), (e._unclipped = e.slice(0)));
            for (var n = 0; n <= 3; n++)
              n < 3
                ? ((e[n] < 0 || e[n] > 255) && (e._clipped = !0),
                  (e[n] = t(e[n], 0, 255)))
                : n === 3 && (e[n] = t(e[n], 0, 1));
            return e;
          },
          r = {},
          i = 0,
          a = [
            `Boolean`,
            `Number`,
            `String`,
            `Function`,
            `Array`,
            `Date`,
            `RegExp`,
            `Undefined`,
            `Null`,
          ];
        i < a.length;
        i += 1
      ) {
        var o = a[i];
        r[`[object ` + o + `]`] = o.toLowerCase();
      }
      var s = function (e) {
          return r[Object.prototype.toString.call(e)] || `object`;
        },
        c = s,
        l = function (e, t) {
          return (
            t === void 0 && (t = null),
            e.length >= 3
              ? Array.prototype.slice.call(e)
              : c(e[0]) == `object` && t
                ? t
                    .split(``)
                    .filter(function (t) {
                      return e[0][t] !== void 0;
                    })
                    .map(function (t) {
                      return e[0][t];
                    })
                : e[0]
          );
        },
        u = s,
        d = function (e) {
          if (e.length < 2) return null;
          var t = e.length - 1;
          return u(e[t]) == `string` ? e[t].toLowerCase() : null;
        },
        f = Math.PI,
        p = {
          clip_rgb: n,
          limit: e,
          type: s,
          unpack: l,
          last: d,
          PI: f,
          TWOPI: f * 2,
          PITHIRD: f / 3,
          DEG2RAD: f / 180,
          RAD2DEG: 180 / f,
        },
        m = { format: {}, autodetect: [] },
        h = p.last,
        g = p.clip_rgb,
        _ = p.type,
        v = m,
        y = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          var n = this;
          if (
            _(e[0]) === `object` &&
            e[0].constructor &&
            e[0].constructor === this.constructor
          )
            return e[0];
          var r = h(e),
            i = !1;
          if (!r) {
            ((i = !0),
              (v.sorted ||=
                ((v.autodetect = v.autodetect.sort(function (e, t) {
                  return t.p - e.p;
                })),
                !0)));
            for (var a = 0, o = v.autodetect; a < o.length; a += 1) {
              var s = o[a];
              if (((r = s.test.apply(s, e)), r)) break;
            }
          }
          if (v.format[r])
            n._rgb = g(v.format[r].apply(null, i ? e : e.slice(0, -1)));
          else throw Error(`unknown format: ` + e);
          n._rgb.length === 3 && n._rgb.push(1);
        };
      y.prototype.toString = function () {
        return _(this.hex) == `function`
          ? this.hex()
          : `[` + this._rgb.join(`,`) + `]`;
      };
      var b = y,
        x = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          return new (Function.prototype.bind.apply(
            x.Color,
            [null].concat(e),
          ))();
        };
      ((x.Color = b), (x.version = `2.4.2`));
      var S = x,
        C = p.unpack,
        w = Math.max,
        T = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          var n = C(e, `rgb`),
            r = n[0],
            i = n[1],
            a = n[2];
          ((r /= 255), (i /= 255), (a /= 255));
          var o = 1 - w(r, w(i, a)),
            s = o < 1 ? 1 / (1 - o) : 0;
          return [(1 - r - o) * s, (1 - i - o) * s, (1 - a - o) * s, o];
        },
        E = p.unpack,
        D = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          e = E(e, `cmyk`);
          var n = e[0],
            r = e[1],
            i = e[2],
            a = e[3],
            o = e.length > 4 ? e[4] : 1;
          return a === 1
            ? [0, 0, 0, o]
            : [
                n >= 1 ? 0 : 255 * (1 - n) * (1 - a),
                r >= 1 ? 0 : 255 * (1 - r) * (1 - a),
                i >= 1 ? 0 : 255 * (1 - i) * (1 - a),
                o,
              ];
        },
        O = S,
        k = b,
        A = m,
        ee = p.unpack,
        j = p.type,
        te = T;
      ((k.prototype.cmyk = function () {
        return te(this._rgb);
      }),
        (O.cmyk = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          return new (Function.prototype.bind.apply(
            k,
            [null].concat(e, [`cmyk`]),
          ))();
        }),
        (A.format.cmyk = D),
        A.autodetect.push({
          p: 2,
          test: function () {
            for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
            if (((e = ee(e, `cmyk`)), j(e) === `array` && e.length === 4))
              return `cmyk`;
          },
        }));
      var M = p.unpack,
        N = p.last,
        P = function (e) {
          return Math.round(e * 100) / 100;
        },
        ne = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          var n = M(e, `hsla`),
            r = N(e) || `lsa`;
          return (
            (n[0] = P(n[0] || 0)),
            (n[1] = P(n[1] * 100) + `%`),
            (n[2] = P(n[2] * 100) + `%`),
            r === `hsla` || (n.length > 3 && n[3] < 1)
              ? ((n[3] = n.length > 3 ? n[3] : 1), (r = `hsla`))
              : (n.length = 3),
            r + `(` + n.join(`,`) + `)`
          );
        },
        F = p.unpack,
        I = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          e = F(e, `rgba`);
          var n = e[0],
            r = e[1],
            i = e[2];
          ((n /= 255), (r /= 255), (i /= 255));
          var a = Math.min(n, r, i),
            o = Math.max(n, r, i),
            s = (o + a) / 2,
            c,
            l;
          return (
            o === a
              ? ((c = 0), (l = NaN))
              : (c = s < 0.5 ? (o - a) / (o + a) : (o - a) / (2 - o - a)),
            n == o
              ? (l = (r - i) / (o - a))
              : r == o
                ? (l = 2 + (i - n) / (o - a))
                : i == o && (l = 4 + (n - r) / (o - a)),
            (l *= 60),
            l < 0 && (l += 360),
            e.length > 3 && e[3] !== void 0 ? [l, c, s, e[3]] : [l, c, s]
          );
        },
        re = p.unpack,
        ie = p.last,
        ae = ne,
        oe = I,
        se = Math.round,
        ce = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          var n = re(e, `rgba`),
            r = ie(e) || `rgb`;
          return r.substr(0, 3) == `hsl`
            ? ae(oe(n), r)
            : ((n[0] = se(n[0])),
              (n[1] = se(n[1])),
              (n[2] = se(n[2])),
              (r === `rgba` || (n.length > 3 && n[3] < 1)) &&
                ((n[3] = n.length > 3 ? n[3] : 1), (r = `rgba`)),
              r + `(` + n.slice(0, r === `rgb` ? 3 : 4).join(`,`) + `)`);
        },
        le = p.unpack,
        ue = Math.round,
        L = function () {
          for (var e, t = [], n = arguments.length; n--;) t[n] = arguments[n];
          t = le(t, `hsl`);
          var r = t[0],
            i = t[1],
            a = t[2],
            o,
            s,
            c;
          if (i === 0) o = s = c = a * 255;
          else {
            var l = [0, 0, 0],
              u = [0, 0, 0],
              d = a < 0.5 ? a * (1 + i) : a + i - a * i,
              f = 2 * a - d,
              p = r / 360;
            ((l[0] = p + 1 / 3), (l[1] = p), (l[2] = p - 1 / 3));
            for (var m = 0; m < 3; m++)
              (l[m] < 0 && (l[m] += 1),
                l[m] > 1 && --l[m],
                6 * l[m] < 1
                  ? (u[m] = f + (d - f) * 6 * l[m])
                  : 2 * l[m] < 1
                    ? (u[m] = d)
                    : 3 * l[m] < 2
                      ? (u[m] = f + (d - f) * (2 / 3 - l[m]) * 6)
                      : (u[m] = f));
            ((e = [ue(u[0] * 255), ue(u[1] * 255), ue(u[2] * 255)]),
              (o = e[0]),
              (s = e[1]),
              (c = e[2]));
          }
          return t.length > 3 ? [o, s, c, t[3]] : [o, s, c, 1];
        },
        R = L,
        z = m,
        B = /^rgb\(\s*(-?\d+),\s*(-?\d+)\s*,\s*(-?\d+)\s*\)$/,
        de =
          /^rgba\(\s*(-?\d+),\s*(-?\d+)\s*,\s*(-?\d+)\s*,\s*([01]|[01]?\.\d+)\)$/,
        fe =
          /^rgb\(\s*(-?\d+(?:\.\d+)?)%,\s*(-?\d+(?:\.\d+)?)%\s*,\s*(-?\d+(?:\.\d+)?)%\s*\)$/,
        pe =
          /^rgba\(\s*(-?\d+(?:\.\d+)?)%,\s*(-?\d+(?:\.\d+)?)%\s*,\s*(-?\d+(?:\.\d+)?)%\s*,\s*([01]|[01]?\.\d+)\)$/,
        me =
          /^hsl\(\s*(-?\d+(?:\.\d+)?),\s*(-?\d+(?:\.\d+)?)%\s*,\s*(-?\d+(?:\.\d+)?)%\s*\)$/,
        he =
          /^hsla\(\s*(-?\d+(?:\.\d+)?),\s*(-?\d+(?:\.\d+)?)%\s*,\s*(-?\d+(?:\.\d+)?)%\s*,\s*([01]|[01]?\.\d+)\)$/,
        ge = Math.round,
        V = function (e) {
          e = e.toLowerCase().trim();
          var t;
          if (z.format.named)
            try {
              return z.format.named(e);
            } catch {}
          if ((t = e.match(B))) {
            for (var n = t.slice(1, 4), r = 0; r < 3; r++) n[r] = +n[r];
            return ((n[3] = 1), n);
          }
          if ((t = e.match(de))) {
            for (var i = t.slice(1, 5), a = 0; a < 4; a++) i[a] = +i[a];
            return i;
          }
          if ((t = e.match(fe))) {
            for (var o = t.slice(1, 4), s = 0; s < 3; s++)
              o[s] = ge(o[s] * 2.55);
            return ((o[3] = 1), o);
          }
          if ((t = e.match(pe))) {
            for (var c = t.slice(1, 5), l = 0; l < 3; l++)
              c[l] = ge(c[l] * 2.55);
            return ((c[3] = +c[3]), c);
          }
          if ((t = e.match(me))) {
            var u = t.slice(1, 4);
            ((u[1] *= 0.01), (u[2] *= 0.01));
            var d = R(u);
            return ((d[3] = 1), d);
          }
          if ((t = e.match(he))) {
            var f = t.slice(1, 4);
            ((f[1] *= 0.01), (f[2] *= 0.01));
            var p = R(f);
            return ((p[3] = +t[4]), p);
          }
        };
      V.test = function (e) {
        return (
          B.test(e) ||
          de.test(e) ||
          fe.test(e) ||
          pe.test(e) ||
          me.test(e) ||
          he.test(e)
        );
      };
      var _e = V,
        ve = S,
        ye = b,
        be = m,
        xe = p.type,
        Se = ce,
        Ce = _e;
      ((ye.prototype.css = function (e) {
        return Se(this._rgb, e);
      }),
        (ve.css = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          return new (Function.prototype.bind.apply(
            ye,
            [null].concat(e, [`css`]),
          ))();
        }),
        (be.format.css = Ce),
        be.autodetect.push({
          p: 5,
          test: function (e) {
            for (var t = [], n = arguments.length - 1; n-- > 0;)
              t[n] = arguments[n + 1];
            if (!t.length && xe(e) === `string` && Ce.test(e)) return `css`;
          },
        }));
      var we = b,
        Te = S,
        Ee = m,
        De = p.unpack;
      ((Ee.format.gl = function () {
        for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
        var n = De(e, `rgba`);
        return ((n[0] *= 255), (n[1] *= 255), (n[2] *= 255), n);
      }),
        (Te.gl = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          return new (Function.prototype.bind.apply(
            we,
            [null].concat(e, [`gl`]),
          ))();
        }),
        (we.prototype.gl = function () {
          var e = this._rgb;
          return [e[0] / 255, e[1] / 255, e[2] / 255, e[3]];
        }));
      var Oe = p.unpack,
        ke = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          var n = Oe(e, `rgb`),
            r = n[0],
            i = n[1],
            a = n[2],
            o = Math.min(r, i, a),
            s = Math.max(r, i, a),
            c = s - o,
            l = (c * 100) / 255,
            u = (o / (255 - c)) * 100,
            d;
          return (
            c === 0
              ? (d = NaN)
              : (r === s && (d = (i - a) / c),
                i === s && (d = 2 + (a - r) / c),
                a === s && (d = 4 + (r - i) / c),
                (d *= 60),
                d < 0 && (d += 360)),
            [d, l, u]
          );
        },
        Ae = p.unpack,
        je = Math.floor,
        Me = function () {
          for (var e, t, n, r, i, a, o = [], s = arguments.length; s--;)
            o[s] = arguments[s];
          o = Ae(o, `hcg`);
          var c = o[0],
            l = o[1],
            u = o[2],
            d,
            f,
            p;
          u *= 255;
          var m = l * 255;
          if (l === 0) d = f = p = u;
          else {
            (c === 360 && (c = 0),
              c > 360 && (c -= 360),
              c < 0 && (c += 360),
              (c /= 60));
            var h = je(c),
              g = c - h,
              _ = u * (1 - l),
              v = _ + m * (1 - g),
              y = _ + m * g,
              b = _ + m;
            switch (h) {
              case 0:
                ((e = [b, y, _]), (d = e[0]), (f = e[1]), (p = e[2]));
                break;
              case 1:
                ((t = [v, b, _]), (d = t[0]), (f = t[1]), (p = t[2]));
                break;
              case 2:
                ((n = [_, b, y]), (d = n[0]), (f = n[1]), (p = n[2]));
                break;
              case 3:
                ((r = [_, v, b]), (d = r[0]), (f = r[1]), (p = r[2]));
                break;
              case 4:
                ((i = [y, _, b]), (d = i[0]), (f = i[1]), (p = i[2]));
                break;
              case 5:
                ((a = [b, _, v]), (d = a[0]), (f = a[1]), (p = a[2]));
            }
          }
          return [d, f, p, o.length > 3 ? o[3] : 1];
        },
        Ne = p.unpack,
        Pe = p.type,
        Fe = S,
        Ie = b,
        Le = m,
        Re = ke;
      ((Ie.prototype.hcg = function () {
        return Re(this._rgb);
      }),
        (Fe.hcg = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          return new (Function.prototype.bind.apply(
            Ie,
            [null].concat(e, [`hcg`]),
          ))();
        }),
        (Le.format.hcg = Me),
        Le.autodetect.push({
          p: 1,
          test: function () {
            for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
            if (((e = Ne(e, `hcg`)), Pe(e) === `array` && e.length === 3))
              return `hcg`;
          },
        }));
      var ze = p.unpack,
        Be = p.last,
        Ve = Math.round,
        He = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          var n = ze(e, `rgba`),
            r = n[0],
            i = n[1],
            a = n[2],
            o = n[3],
            s = Be(e) || `auto`;
          (o === void 0 && (o = 1),
            s === `auto` && (s = o < 1 ? `rgba` : `rgb`),
            (r = Ve(r)),
            (i = Ve(i)),
            (a = Ve(a)));
          var c = `000000` + ((r << 16) | (i << 8) | a).toString(16);
          c = c.substr(c.length - 6);
          var l = `0` + Ve(o * 255).toString(16);
          switch (((l = l.substr(l.length - 2)), s.toLowerCase())) {
            case `rgba`:
              return `#` + c + l;
            case `argb`:
              return `#` + l + c;
            default:
              return `#` + c;
          }
        },
        Ue = /^#?([A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$/,
        H = /^#?([A-Fa-f0-9]{8}|[A-Fa-f0-9]{4})$/,
        We = function (e) {
          if (e.match(Ue)) {
            ((e.length === 4 || e.length === 7) && (e = e.substr(1)),
              e.length === 3 &&
                ((e = e.split(``)),
                (e = e[0] + e[0] + e[1] + e[1] + e[2] + e[2])));
            var t = parseInt(e, 16);
            return [t >> 16, (t >> 8) & 255, t & 255, 1];
          }
          if (e.match(H)) {
            ((e.length === 5 || e.length === 9) && (e = e.substr(1)),
              e.length === 4 &&
                ((e = e.split(``)),
                (e = e[0] + e[0] + e[1] + e[1] + e[2] + e[2] + e[3] + e[3])));
            var n = parseInt(e, 16);
            return [
              (n >> 24) & 255,
              (n >> 16) & 255,
              (n >> 8) & 255,
              Math.round(((n & 255) / 255) * 100) / 100,
            ];
          }
          throw Error(`unknown hex color: ` + e);
        },
        U = S,
        Ge = b,
        W = p.type,
        Ke = m,
        G = He;
      ((Ge.prototype.hex = function (e) {
        return G(this._rgb, e);
      }),
        (U.hex = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          return new (Function.prototype.bind.apply(
            Ge,
            [null].concat(e, [`hex`]),
          ))();
        }),
        (Ke.format.hex = We),
        Ke.autodetect.push({
          p: 4,
          test: function (e) {
            for (var t = [], n = arguments.length - 1; n-- > 0;)
              t[n] = arguments[n + 1];
            if (
              !t.length &&
              W(e) === `string` &&
              [3, 4, 5, 6, 7, 8, 9].indexOf(e.length) >= 0
            )
              return `hex`;
          },
        }));
      var K = p.unpack,
        qe = p.TWOPI,
        Je = Math.min,
        Ye = Math.sqrt,
        q = Math.acos,
        Xe = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          var n = K(e, `rgb`),
            r = n[0],
            i = n[1],
            a = n[2];
          ((r /= 255), (i /= 255), (a /= 255));
          var o,
            s = Je(r, i, a),
            c = (r + i + a) / 3,
            l = c > 0 ? 1 - s / c : 0;
          return (
            l === 0
              ? (o = NaN)
              : ((o = (r - i + (r - a)) / 2),
                (o /= Ye((r - i) * (r - i) + (r - a) * (i - a))),
                (o = q(o)),
                a > i && (o = qe - o),
                (o /= qe)),
            [o * 360, l, c]
          );
        },
        Ze = p.unpack,
        Qe = p.limit,
        J = p.TWOPI,
        $e = p.PITHIRD,
        et = Math.cos,
        tt = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          e = Ze(e, `hsi`);
          var n = e[0],
            r = e[1],
            i = e[2],
            a,
            o,
            s;
          return (
            isNaN(n) && (n = 0),
            isNaN(r) && (r = 0),
            n > 360 && (n -= 360),
            n < 0 && (n += 360),
            (n /= 360),
            n < 1 / 3
              ? ((s = (1 - r) / 3),
                (a = (1 + (r * et(J * n)) / et($e - J * n)) / 3),
                (o = 1 - (s + a)))
              : n < 2 / 3
                ? ((n -= 1 / 3),
                  (a = (1 - r) / 3),
                  (o = (1 + (r * et(J * n)) / et($e - J * n)) / 3),
                  (s = 1 - (a + o)))
                : ((n -= 2 / 3),
                  (o = (1 - r) / 3),
                  (s = (1 + (r * et(J * n)) / et($e - J * n)) / 3),
                  (a = 1 - (o + s))),
            (a = Qe(i * a * 3)),
            (o = Qe(i * o * 3)),
            (s = Qe(i * s * 3)),
            [a * 255, o * 255, s * 255, e.length > 3 ? e[3] : 1]
          );
        },
        nt = p.unpack,
        rt = p.type,
        it = S,
        at = b,
        ot = m,
        st = Xe;
      ((at.prototype.hsi = function () {
        return st(this._rgb);
      }),
        (it.hsi = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          return new (Function.prototype.bind.apply(
            at,
            [null].concat(e, [`hsi`]),
          ))();
        }),
        (ot.format.hsi = tt),
        ot.autodetect.push({
          p: 2,
          test: function () {
            for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
            if (((e = nt(e, `hsi`)), rt(e) === `array` && e.length === 3))
              return `hsi`;
          },
        }));
      var ct = p.unpack,
        lt = p.type,
        ut = S,
        dt = b,
        ft = m,
        pt = I;
      ((dt.prototype.hsl = function () {
        return pt(this._rgb);
      }),
        (ut.hsl = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          return new (Function.prototype.bind.apply(
            dt,
            [null].concat(e, [`hsl`]),
          ))();
        }),
        (ft.format.hsl = L),
        ft.autodetect.push({
          p: 2,
          test: function () {
            for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
            if (((e = ct(e, `hsl`)), lt(e) === `array` && e.length === 3))
              return `hsl`;
          },
        }));
      var mt = p.unpack,
        ht = Math.min,
        gt = Math.max,
        _t = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          e = mt(e, `rgb`);
          var n = e[0],
            r = e[1],
            i = e[2],
            a = ht(n, r, i),
            o = gt(n, r, i),
            s = o - a,
            c,
            l,
            u;
          return (
            (u = o / 255),
            o === 0
              ? ((c = NaN), (l = 0))
              : ((l = s / o),
                n === o && (c = (r - i) / s),
                r === o && (c = 2 + (i - n) / s),
                i === o && (c = 4 + (n - r) / s),
                (c *= 60),
                c < 0 && (c += 360)),
            [c, l, u]
          );
        },
        vt = p.unpack,
        yt = Math.floor,
        bt = function () {
          for (var e, t, n, r, i, a, o = [], s = arguments.length; s--;)
            o[s] = arguments[s];
          o = vt(o, `hsv`);
          var c = o[0],
            l = o[1],
            u = o[2],
            d,
            f,
            p;
          if (((u *= 255), l === 0)) d = f = p = u;
          else {
            (c === 360 && (c = 0),
              c > 360 && (c -= 360),
              c < 0 && (c += 360),
              (c /= 60));
            var m = yt(c),
              h = c - m,
              g = u * (1 - l),
              _ = u * (1 - l * h),
              v = u * (1 - l * (1 - h));
            switch (m) {
              case 0:
                ((e = [u, v, g]), (d = e[0]), (f = e[1]), (p = e[2]));
                break;
              case 1:
                ((t = [_, u, g]), (d = t[0]), (f = t[1]), (p = t[2]));
                break;
              case 2:
                ((n = [g, u, v]), (d = n[0]), (f = n[1]), (p = n[2]));
                break;
              case 3:
                ((r = [g, _, u]), (d = r[0]), (f = r[1]), (p = r[2]));
                break;
              case 4:
                ((i = [v, g, u]), (d = i[0]), (f = i[1]), (p = i[2]));
                break;
              case 5:
                ((a = [u, g, _]), (d = a[0]), (f = a[1]), (p = a[2]));
            }
          }
          return [d, f, p, o.length > 3 ? o[3] : 1];
        },
        xt = p.unpack,
        St = p.type,
        Ct = S,
        wt = b,
        Tt = m,
        Et = _t;
      ((wt.prototype.hsv = function () {
        return Et(this._rgb);
      }),
        (Ct.hsv = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          return new (Function.prototype.bind.apply(
            wt,
            [null].concat(e, [`hsv`]),
          ))();
        }),
        (Tt.format.hsv = bt),
        Tt.autodetect.push({
          p: 2,
          test: function () {
            for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
            if (((e = xt(e, `hsv`)), St(e) === `array` && e.length === 3))
              return `hsv`;
          },
        }));
      var Dt = {
          Kn: 18,
          Xn: 0.95047,
          Yn: 1,
          Zn: 1.08883,
          t0: 0.137931034,
          t1: 0.206896552,
          t2: 0.12841855,
          t3: 0.008856452,
        },
        Ot = Dt,
        kt = p.unpack,
        At = Math.pow,
        jt = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          var n = kt(e, `rgb`),
            r = n[0],
            i = n[1],
            a = n[2],
            o = Pt(r, i, a),
            s = o[0],
            c = o[1],
            l = o[2],
            u = 116 * c - 16;
          return [u < 0 ? 0 : u, 500 * (s - c), 200 * (c - l)];
        },
        Mt = function (e) {
          return (e /= 255) <= 0.04045
            ? e / 12.92
            : At((e + 0.055) / 1.055, 2.4);
        },
        Nt = function (e) {
          return e > Ot.t3 ? At(e, 1 / 3) : e / Ot.t2 + Ot.t0;
        },
        Pt = function (e, t, n) {
          return (
            (e = Mt(e)),
            (t = Mt(t)),
            (n = Mt(n)),
            [
              Nt((0.4124564 * e + 0.3575761 * t + 0.1804375 * n) / Ot.Xn),
              Nt((0.2126729 * e + 0.7151522 * t + 0.072175 * n) / Ot.Yn),
              Nt((0.0193339 * e + 0.119192 * t + 0.9503041 * n) / Ot.Zn),
            ]
          );
        },
        Ft = jt,
        It = Dt,
        Lt = p.unpack,
        Rt = Math.pow,
        zt = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          e = Lt(e, `lab`);
          var n = e[0],
            r = e[1],
            i = e[2],
            a,
            o,
            s,
            c,
            l,
            u;
          return (
            (o = (n + 16) / 116),
            (a = isNaN(r) ? o : o + r / 500),
            (s = isNaN(i) ? o : o - i / 200),
            (o = It.Yn * Vt(o)),
            (a = It.Xn * Vt(a)),
            (s = It.Zn * Vt(s)),
            (c = Bt(3.2404542 * a - 1.5371385 * o - 0.4985314 * s)),
            (l = Bt(-0.969266 * a + 1.8760108 * o + 0.041556 * s)),
            (u = Bt(0.0556434 * a - 0.2040259 * o + 1.0572252 * s)),
            [c, l, u, e.length > 3 ? e[3] : 1]
          );
        },
        Bt = function (e) {
          return (
            255 * (e <= 0.00304 ? 12.92 * e : 1.055 * Rt(e, 1 / 2.4) - 0.055)
          );
        },
        Vt = function (e) {
          return e > It.t1 ? e * e * e : It.t2 * (e - It.t0);
        },
        Ht = zt,
        Ut = p.unpack,
        Wt = p.type,
        Gt = S,
        Kt = b,
        qt = m,
        Jt = Ft;
      ((Kt.prototype.lab = function () {
        return Jt(this._rgb);
      }),
        (Gt.lab = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          return new (Function.prototype.bind.apply(
            Kt,
            [null].concat(e, [`lab`]),
          ))();
        }),
        (qt.format.lab = Ht),
        qt.autodetect.push({
          p: 2,
          test: function () {
            for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
            if (((e = Ut(e, `lab`)), Wt(e) === `array` && e.length === 3))
              return `lab`;
          },
        }));
      var Yt = p.unpack,
        Xt = p.RAD2DEG,
        Zt = Math.sqrt,
        Qt = Math.atan2,
        $t = Math.round,
        en = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          var n = Yt(e, `lab`),
            r = n[0],
            i = n[1],
            a = n[2],
            o = Zt(i * i + a * a),
            s = (Qt(a, i) * Xt + 360) % 360;
          return ($t(o * 1e4) === 0 && (s = NaN), [r, o, s]);
        },
        tn = p.unpack,
        nn = Ft,
        rn = en,
        an = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          var n = tn(e, `rgb`),
            r = n[0],
            i = n[1],
            a = n[2],
            o = nn(r, i, a),
            s = o[0],
            c = o[1],
            l = o[2];
          return rn(s, c, l);
        },
        on = p.unpack,
        sn = p.DEG2RAD,
        cn = Math.sin,
        ln = Math.cos,
        un = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          var n = on(e, `lch`),
            r = n[0],
            i = n[1],
            a = n[2];
          return (isNaN(a) && (a = 0), (a *= sn), [r, ln(a) * i, cn(a) * i]);
        },
        dn = p.unpack,
        fn = un,
        pn = Ht,
        mn = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          e = dn(e, `lch`);
          var n = e[0],
            r = e[1],
            i = e[2],
            a = fn(n, r, i),
            o = a[0],
            s = a[1],
            c = a[2],
            l = pn(o, s, c);
          return [l[0], l[1], l[2], e.length > 3 ? e[3] : 1];
        },
        hn = p.unpack,
        gn = mn,
        _n = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          var n = hn(e, `hcl`).reverse();
          return gn.apply(void 0, n);
        },
        vn = p.unpack,
        yn = p.type,
        bn = S,
        xn = b,
        Sn = m,
        Cn = an;
      ((xn.prototype.lch = function () {
        return Cn(this._rgb);
      }),
        (xn.prototype.hcl = function () {
          return Cn(this._rgb).reverse();
        }),
        (bn.lch = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          return new (Function.prototype.bind.apply(
            xn,
            [null].concat(e, [`lch`]),
          ))();
        }),
        (bn.hcl = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          return new (Function.prototype.bind.apply(
            xn,
            [null].concat(e, [`hcl`]),
          ))();
        }),
        (Sn.format.lch = mn),
        (Sn.format.hcl = _n),
        [`lch`, `hcl`].forEach(function (e) {
          return Sn.autodetect.push({
            p: 2,
            test: function () {
              for (var t = [], n = arguments.length; n--;) t[n] = arguments[n];
              if (((t = vn(t, e)), yn(t) === `array` && t.length === 3))
                return e;
            },
          });
        }));
      var wn = {
          aliceblue: `#f0f8ff`,
          antiquewhite: `#faebd7`,
          aqua: `#00ffff`,
          aquamarine: `#7fffd4`,
          azure: `#f0ffff`,
          beige: `#f5f5dc`,
          bisque: `#ffe4c4`,
          black: `#000000`,
          blanchedalmond: `#ffebcd`,
          blue: `#0000ff`,
          blueviolet: `#8a2be2`,
          brown: `#a52a2a`,
          burlywood: `#deb887`,
          cadetblue: `#5f9ea0`,
          chartreuse: `#7fff00`,
          chocolate: `#d2691e`,
          coral: `#ff7f50`,
          cornflower: `#6495ed`,
          cornflowerblue: `#6495ed`,
          cornsilk: `#fff8dc`,
          crimson: `#dc143c`,
          cyan: `#00ffff`,
          darkblue: `#00008b`,
          darkcyan: `#008b8b`,
          darkgoldenrod: `#b8860b`,
          darkgray: `#a9a9a9`,
          darkgreen: `#006400`,
          darkgrey: `#a9a9a9`,
          darkkhaki: `#bdb76b`,
          darkmagenta: `#8b008b`,
          darkolivegreen: `#556b2f`,
          darkorange: `#ff8c00`,
          darkorchid: `#9932cc`,
          darkred: `#8b0000`,
          darksalmon: `#e9967a`,
          darkseagreen: `#8fbc8f`,
          darkslateblue: `#483d8b`,
          darkslategray: `#2f4f4f`,
          darkslategrey: `#2f4f4f`,
          darkturquoise: `#00ced1`,
          darkviolet: `#9400d3`,
          deeppink: `#ff1493`,
          deepskyblue: `#00bfff`,
          dimgray: `#696969`,
          dimgrey: `#696969`,
          dodgerblue: `#1e90ff`,
          firebrick: `#b22222`,
          floralwhite: `#fffaf0`,
          forestgreen: `#228b22`,
          fuchsia: `#ff00ff`,
          gainsboro: `#dcdcdc`,
          ghostwhite: `#f8f8ff`,
          gold: `#ffd700`,
          goldenrod: `#daa520`,
          gray: `#808080`,
          green: `#008000`,
          greenyellow: `#adff2f`,
          grey: `#808080`,
          honeydew: `#f0fff0`,
          hotpink: `#ff69b4`,
          indianred: `#cd5c5c`,
          indigo: `#4b0082`,
          ivory: `#fffff0`,
          khaki: `#f0e68c`,
          laserlemon: `#ffff54`,
          lavender: `#e6e6fa`,
          lavenderblush: `#fff0f5`,
          lawngreen: `#7cfc00`,
          lemonchiffon: `#fffacd`,
          lightblue: `#add8e6`,
          lightcoral: `#f08080`,
          lightcyan: `#e0ffff`,
          lightgoldenrod: `#fafad2`,
          lightgoldenrodyellow: `#fafad2`,
          lightgray: `#d3d3d3`,
          lightgreen: `#90ee90`,
          lightgrey: `#d3d3d3`,
          lightpink: `#ffb6c1`,
          lightsalmon: `#ffa07a`,
          lightseagreen: `#20b2aa`,
          lightskyblue: `#87cefa`,
          lightslategray: `#778899`,
          lightslategrey: `#778899`,
          lightsteelblue: `#b0c4de`,
          lightyellow: `#ffffe0`,
          lime: `#00ff00`,
          limegreen: `#32cd32`,
          linen: `#faf0e6`,
          magenta: `#ff00ff`,
          maroon: `#800000`,
          maroon2: `#7f0000`,
          maroon3: `#b03060`,
          mediumaquamarine: `#66cdaa`,
          mediumblue: `#0000cd`,
          mediumorchid: `#ba55d3`,
          mediumpurple: `#9370db`,
          mediumseagreen: `#3cb371`,
          mediumslateblue: `#7b68ee`,
          mediumspringgreen: `#00fa9a`,
          mediumturquoise: `#48d1cc`,
          mediumvioletred: `#c71585`,
          midnightblue: `#191970`,
          mintcream: `#f5fffa`,
          mistyrose: `#ffe4e1`,
          moccasin: `#ffe4b5`,
          navajowhite: `#ffdead`,
          navy: `#000080`,
          oldlace: `#fdf5e6`,
          olive: `#808000`,
          olivedrab: `#6b8e23`,
          orange: `#ffa500`,
          orangered: `#ff4500`,
          orchid: `#da70d6`,
          palegoldenrod: `#eee8aa`,
          palegreen: `#98fb98`,
          paleturquoise: `#afeeee`,
          palevioletred: `#db7093`,
          papayawhip: `#ffefd5`,
          peachpuff: `#ffdab9`,
          peru: `#cd853f`,
          pink: `#ffc0cb`,
          plum: `#dda0dd`,
          powderblue: `#b0e0e6`,
          purple: `#800080`,
          purple2: `#7f007f`,
          purple3: `#a020f0`,
          rebeccapurple: `#663399`,
          red: `#ff0000`,
          rosybrown: `#bc8f8f`,
          royalblue: `#4169e1`,
          saddlebrown: `#8b4513`,
          salmon: `#fa8072`,
          sandybrown: `#f4a460`,
          seagreen: `#2e8b57`,
          seashell: `#fff5ee`,
          sienna: `#a0522d`,
          silver: `#c0c0c0`,
          skyblue: `#87ceeb`,
          slateblue: `#6a5acd`,
          slategray: `#708090`,
          slategrey: `#708090`,
          snow: `#fffafa`,
          springgreen: `#00ff7f`,
          steelblue: `#4682b4`,
          tan: `#d2b48c`,
          teal: `#008080`,
          thistle: `#d8bfd8`,
          tomato: `#ff6347`,
          turquoise: `#40e0d0`,
          violet: `#ee82ee`,
          wheat: `#f5deb3`,
          white: `#ffffff`,
          whitesmoke: `#f5f5f5`,
          yellow: `#ffff00`,
          yellowgreen: `#9acd32`,
        },
        Tn = b,
        En = m,
        Dn = p.type,
        On = wn,
        kn = We,
        An = He;
      ((Tn.prototype.name = function () {
        for (
          var e = An(this._rgb, `rgb`), t = 0, n = Object.keys(On);
          t < n.length;
          t += 1
        ) {
          var r = n[t];
          if (On[r] === e) return r.toLowerCase();
        }
        return e;
      }),
        (En.format.named = function (e) {
          if (((e = e.toLowerCase()), On[e])) return kn(On[e]);
          throw Error(`unknown color name: ` + e);
        }),
        En.autodetect.push({
          p: 5,
          test: function (e) {
            for (var t = [], n = arguments.length - 1; n-- > 0;)
              t[n] = arguments[n + 1];
            if (!t.length && Dn(e) === `string` && On[e.toLowerCase()])
              return `named`;
          },
        }));
      var jn = p.unpack,
        Mn = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          var n = jn(e, `rgb`),
            r = n[0],
            i = n[1],
            a = n[2];
          return (r << 16) + (i << 8) + a;
        },
        Nn = p.type,
        Pn = function (e) {
          if (Nn(e) == `number` && e >= 0 && e <= 16777215)
            return [e >> 16, (e >> 8) & 255, e & 255, 1];
          throw Error(`unknown num color: ` + e);
        },
        Fn = S,
        In = b,
        Ln = m,
        Rn = p.type,
        zn = Mn;
      ((In.prototype.num = function () {
        return zn(this._rgb);
      }),
        (Fn.num = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          return new (Function.prototype.bind.apply(
            In,
            [null].concat(e, [`num`]),
          ))();
        }),
        (Ln.format.num = Pn),
        Ln.autodetect.push({
          p: 5,
          test: function () {
            for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
            if (
              e.length === 1 &&
              Rn(e[0]) === `number` &&
              e[0] >= 0 &&
              e[0] <= 16777215
            )
              return `num`;
          },
        }));
      var Bn = S,
        Vn = b,
        Hn = m,
        Un = p.unpack,
        Wn = p.type,
        Gn = Math.round;
      ((Vn.prototype.rgb = function (e) {
        return (
          e === void 0 && (e = !0),
          e === !1 ? this._rgb.slice(0, 3) : this._rgb.slice(0, 3).map(Gn)
        );
      }),
        (Vn.prototype.rgba = function (e) {
          return (
            e === void 0 && (e = !0),
            this._rgb.slice(0, 4).map(function (t, n) {
              return n < 3 ? (e === !1 ? t : Gn(t)) : t;
            })
          );
        }),
        (Bn.rgb = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          return new (Function.prototype.bind.apply(
            Vn,
            [null].concat(e, [`rgb`]),
          ))();
        }),
        (Hn.format.rgb = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          var n = Un(e, `rgba`);
          return (n[3] === void 0 && (n[3] = 1), n);
        }),
        Hn.autodetect.push({
          p: 3,
          test: function () {
            for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
            if (
              ((e = Un(e, `rgba`)),
              Wn(e) === `array` &&
                (e.length === 3 ||
                  (e.length === 4 &&
                    Wn(e[3]) == `number` &&
                    e[3] >= 0 &&
                    e[3] <= 1)))
            )
              return `rgb`;
          },
        }));
      var Kn = Math.log,
        qn = function (e) {
          var t = e / 100,
            n,
            r,
            i;
          return (
            t < 66
              ? ((n = 255),
                (r =
                  t < 6
                    ? 0
                    : -155.25485562709179 -
                      0.44596950469579133 * (r = t - 2) +
                      104.49216199393888 * Kn(r)),
                (i =
                  t < 20
                    ? 0
                    : -254.76935184120902 +
                      0.8274096064007395 * (i = t - 10) +
                      115.67994401066147 * Kn(i)))
              : ((n =
                  351.97690566805693 +
                  0.114206453784165 * (n = t - 55) -
                  40.25366309332127 * Kn(n)),
                (r =
                  325.4494125711974 +
                  0.07943456536662342 * (r = t - 50) -
                  28.0852963507957 * Kn(r)),
                (i = 255)),
            [n, r, i, 1]
          );
        },
        Jn = qn,
        Yn = p.unpack,
        Xn = Math.round,
        Zn = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          for (
            var n = Yn(e, `rgb`),
              r = n[0],
              i = n[2],
              a = 1e3,
              o = 4e4,
              s = 0.4,
              c;
            o - a > s;
          ) {
            c = (o + a) * 0.5;
            var l = Jn(c);
            l[2] / l[0] >= i / r ? (o = c) : (a = c);
          }
          return Xn(c);
        },
        Qn = S,
        $n = b,
        er = m,
        tr = Zn;
      (($n.prototype.temp =
        $n.prototype.kelvin =
        $n.prototype.temperature =
          function () {
            return tr(this._rgb);
          }),
        (Qn.temp =
          Qn.kelvin =
          Qn.temperature =
            function () {
              for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
              return new (Function.prototype.bind.apply(
                $n,
                [null].concat(e, [`temp`]),
              ))();
            }),
        (er.format.temp = er.format.kelvin = er.format.temperature = qn));
      var nr = p.unpack,
        rr = Math.cbrt,
        ir = Math.pow,
        ar = Math.sign,
        or = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          var n = nr(e, `rgb`),
            r = n[0],
            i = n[1],
            a = n[2],
            o = [sr(r / 255), sr(i / 255), sr(a / 255)],
            s = o[0],
            c = o[1],
            l = o[2],
            u = rr(0.4122214708 * s + 0.5363325363 * c + 0.0514459929 * l),
            d = rr(0.2119034982 * s + 0.6806995451 * c + 0.1073969566 * l),
            f = rr(0.0883024619 * s + 0.2817188376 * c + 0.6299787005 * l);
          return [
            0.2104542553 * u + 0.793617785 * d - 0.0040720468 * f,
            1.9779984951 * u - 2.428592205 * d + 0.4505937099 * f,
            0.0259040371 * u + 0.7827717662 * d - 0.808675766 * f,
          ];
        };
      function sr(e) {
        var t = Math.abs(e);
        return t < 0.04045
          ? e / 12.92
          : (ar(e) || 1) * ir((t + 0.055) / 1.055, 2.4);
      }
      var cr = p.unpack,
        lr = Math.pow,
        ur = Math.sign,
        dr = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          e = cr(e, `lab`);
          var n = e[0],
            r = e[1],
            i = e[2],
            a = lr(n + 0.3963377774 * r + 0.2158037573 * i, 3),
            o = lr(n - 0.1055613458 * r - 0.0638541728 * i, 3),
            s = lr(n - 0.0894841775 * r - 1.291485548 * i, 3);
          return [
            255 * fr(4.0767416621 * a - 3.3077115913 * o + 0.2309699292 * s),
            255 * fr(-1.2684380046 * a + 2.6097574011 * o - 0.3413193965 * s),
            255 * fr(-0.0041960863 * a - 0.7034186147 * o + 1.707614701 * s),
            e.length > 3 ? e[3] : 1,
          ];
        };
      function fr(e) {
        var t = Math.abs(e);
        return t > 0.0031308
          ? (ur(e) || 1) * (1.055 * lr(t, 1 / 2.4) - 0.055)
          : e * 12.92;
      }
      var pr = p.unpack,
        mr = p.type,
        hr = S,
        gr = b,
        _r = m,
        vr = or;
      ((gr.prototype.oklab = function () {
        return vr(this._rgb);
      }),
        (hr.oklab = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          return new (Function.prototype.bind.apply(
            gr,
            [null].concat(e, [`oklab`]),
          ))();
        }),
        (_r.format.oklab = dr),
        _r.autodetect.push({
          p: 3,
          test: function () {
            for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
            if (((e = pr(e, `oklab`)), mr(e) === `array` && e.length === 3))
              return `oklab`;
          },
        }));
      var yr = p.unpack,
        br = or,
        xr = en,
        Sr = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          var n = yr(e, `rgb`),
            r = n[0],
            i = n[1],
            a = n[2],
            o = br(r, i, a),
            s = o[0],
            c = o[1],
            l = o[2];
          return xr(s, c, l);
        },
        Cr = p.unpack,
        wr = un,
        Tr = dr,
        Er = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          e = Cr(e, `lch`);
          var n = e[0],
            r = e[1],
            i = e[2],
            a = wr(n, r, i),
            o = a[0],
            s = a[1],
            c = a[2],
            l = Tr(o, s, c);
          return [l[0], l[1], l[2], e.length > 3 ? e[3] : 1];
        },
        Dr = p.unpack,
        Or = p.type,
        kr = S,
        Ar = b,
        jr = m,
        Mr = Sr;
      ((Ar.prototype.oklch = function () {
        return Mr(this._rgb);
      }),
        (kr.oklch = function () {
          for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
          return new (Function.prototype.bind.apply(
            Ar,
            [null].concat(e, [`oklch`]),
          ))();
        }),
        (jr.format.oklch = Er),
        jr.autodetect.push({
          p: 3,
          test: function () {
            for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
            if (((e = Dr(e, `oklch`)), Or(e) === `array` && e.length === 3))
              return `oklch`;
          },
        }));
      var Nr = b,
        Pr = p.type;
      Nr.prototype.alpha = function (e, t) {
        return (
          t === void 0 && (t = !1),
          e !== void 0 && Pr(e) === `number`
            ? t
              ? ((this._rgb[3] = e), this)
              : new Nr([this._rgb[0], this._rgb[1], this._rgb[2], e], `rgb`)
            : this._rgb[3]
        );
      };
      var Fr = b;
      Fr.prototype.clipped = function () {
        return this._rgb._clipped || !1;
      };
      var Ir = b,
        Lr = Dt;
      ((Ir.prototype.darken = function (e) {
        e === void 0 && (e = 1);
        var t = this,
          n = t.lab();
        return ((n[0] -= Lr.Kn * e), new Ir(n, `lab`).alpha(t.alpha(), !0));
      }),
        (Ir.prototype.brighten = function (e) {
          return (e === void 0 && (e = 1), this.darken(-e));
        }),
        (Ir.prototype.darker = Ir.prototype.darken),
        (Ir.prototype.brighter = Ir.prototype.brighten));
      var Rr = b;
      Rr.prototype.get = function (e) {
        var t = e.split(`.`),
          n = t[0],
          r = t[1],
          i = this[n]();
        if (r) {
          var a = n.indexOf(r) - (n.substr(0, 2) === `ok` ? 2 : 0);
          if (a > -1) return i[a];
          throw Error(`unknown channel ` + r + ` in mode ` + n);
        }
        return i;
      };
      var zr = b,
        Br = p.type,
        Vr = Math.pow,
        Hr = 1e-7,
        Ur = 20;
      zr.prototype.luminance = function (e) {
        if (e !== void 0 && Br(e) === `number`) {
          if (e === 0) return new zr([0, 0, 0, this._rgb[3]], `rgb`);
          if (e === 1) return new zr([255, 255, 255, this._rgb[3]], `rgb`);
          var t = this.luminance(),
            n = `rgb`,
            r = Ur,
            i = function (t, a) {
              var o = t.interpolate(a, 0.5, n),
                s = o.luminance();
              return Math.abs(e - s) < Hr || !r--
                ? o
                : s > e
                  ? i(t, o)
                  : i(o, a);
            };
          return new zr(
            (t > e
              ? i(new zr([0, 0, 0]), this)
              : i(this, new zr([255, 255, 255]))
            )
              .rgb()
              .concat([this._rgb[3]]),
          );
        }
        return Wr.apply(void 0, this._rgb.slice(0, 3));
      };
      var Wr = function (e, t, n) {
          return (
            (e = Gr(e)),
            (t = Gr(t)),
            (n = Gr(n)),
            0.2126 * e + 0.7152 * t + 0.0722 * n
          );
        },
        Gr = function (e) {
          return (
            (e /= 255),
            e <= 0.03928 ? e / 12.92 : Vr((e + 0.055) / 1.055, 2.4)
          );
        },
        Y = {},
        Kr = b,
        qr = p.type,
        Jr = Y,
        Yr = function (e, t, n) {
          n === void 0 && (n = 0.5);
          for (var r = [], i = arguments.length - 3; i-- > 0;)
            r[i] = arguments[i + 3];
          var a = r[0] || `lrgb`;
          if ((!Jr[a] && !r.length && (a = Object.keys(Jr)[0]), !Jr[a]))
            throw Error(`interpolation mode ` + a + ` is not defined`);
          return (
            qr(e) !== `object` && (e = new Kr(e)),
            qr(t) !== `object` && (t = new Kr(t)),
            Jr[a](e, t, n).alpha(e.alpha() + n * (t.alpha() - e.alpha()))
          );
        },
        Xr = b,
        Zr = Yr;
      Xr.prototype.mix = Xr.prototype.interpolate = function (e, t) {
        t === void 0 && (t = 0.5);
        for (var n = [], r = arguments.length - 2; r-- > 0;)
          n[r] = arguments[r + 2];
        return Zr.apply(void 0, [this, e, t].concat(n));
      };
      var Qr = b;
      Qr.prototype.premultiply = function (e) {
        e === void 0 && (e = !1);
        var t = this._rgb,
          n = t[3];
        return e
          ? ((this._rgb = [t[0] * n, t[1] * n, t[2] * n, n]), this)
          : new Qr([t[0] * n, t[1] * n, t[2] * n, n], `rgb`);
      };
      var $r = b,
        ei = Dt;
      (($r.prototype.saturate = function (e) {
        e === void 0 && (e = 1);
        var t = this,
          n = t.lch();
        return (
          (n[1] += ei.Kn * e),
          n[1] < 0 && (n[1] = 0),
          new $r(n, `lch`).alpha(t.alpha(), !0)
        );
      }),
        ($r.prototype.desaturate = function (e) {
          return (e === void 0 && (e = 1), this.saturate(-e));
        }));
      var ti = b,
        ni = p.type;
      ti.prototype.set = function (e, t, n) {
        n === void 0 && (n = !1);
        var r = e.split(`.`),
          i = r[0],
          a = r[1],
          o = this[i]();
        if (a) {
          var s = i.indexOf(a) - (i.substr(0, 2) === `ok` ? 2 : 0);
          if (s > -1) {
            if (ni(t) == `string`)
              switch (t.charAt(0)) {
                case `+`:
                  o[s] += +t;
                  break;
                case `-`:
                  o[s] += +t;
                  break;
                case `*`:
                  o[s] *= +t.substr(1);
                  break;
                case `/`:
                  o[s] /= +t.substr(1);
                  break;
                default:
                  o[s] = +t;
              }
            else if (ni(t) === `number`) o[s] = t;
            else throw Error(`unsupported value for Color.set`);
            var c = new ti(o, i);
            return n ? ((this._rgb = c._rgb), this) : c;
          }
          throw Error(`unknown channel ` + a + ` in mode ` + i);
        }
        return o;
      };
      var ri = b;
      Y.rgb = function (e, t, n) {
        var r = e._rgb,
          i = t._rgb;
        return new ri(
          r[0] + n * (i[0] - r[0]),
          r[1] + n * (i[1] - r[1]),
          r[2] + n * (i[2] - r[2]),
          `rgb`,
        );
      };
      var ii = b,
        ai = Math.sqrt,
        oi = Math.pow;
      Y.lrgb = function (e, t, n) {
        var r = e._rgb,
          i = r[0],
          a = r[1],
          o = r[2],
          s = t._rgb,
          c = s[0],
          l = s[1],
          u = s[2];
        return new ii(
          ai(oi(i, 2) * (1 - n) + oi(c, 2) * n),
          ai(oi(a, 2) * (1 - n) + oi(l, 2) * n),
          ai(oi(o, 2) * (1 - n) + oi(u, 2) * n),
          `rgb`,
        );
      };
      var si = b;
      Y.lab = function (e, t, n) {
        var r = e.lab(),
          i = t.lab();
        return new si(
          r[0] + n * (i[0] - r[0]),
          r[1] + n * (i[1] - r[1]),
          r[2] + n * (i[2] - r[2]),
          `lab`,
        );
      };
      var ci = b,
        li = function (e, t, n, r) {
          var i, a, o, s;
          r === `hsl`
            ? ((o = e.hsl()), (s = t.hsl()))
            : r === `hsv`
              ? ((o = e.hsv()), (s = t.hsv()))
              : r === `hcg`
                ? ((o = e.hcg()), (s = t.hcg()))
                : r === `hsi`
                  ? ((o = e.hsi()), (s = t.hsi()))
                  : r === `lch` || r === `hcl`
                    ? ((r = `hcl`), (o = e.hcl()), (s = t.hcl()))
                    : r === `oklch` &&
                      ((o = e.oklch().reverse()), (s = t.oklch().reverse()));
          var c, l, u, d, f, p;
          (r.substr(0, 1) === `h` || r === `oklch`) &&
            ((i = o),
            (c = i[0]),
            (u = i[1]),
            (f = i[2]),
            (a = s),
            (l = a[0]),
            (d = a[1]),
            (p = a[2]));
          var m, h, g, _;
          return (
            !isNaN(c) && !isNaN(l)
              ? ((_ =
                  l > c && l - c > 180
                    ? l - (c + 360)
                    : l < c && c - l > 180
                      ? l + 360 - c
                      : l - c),
                (h = c + n * _))
              : isNaN(c)
                ? isNaN(l)
                  ? (h = NaN)
                  : ((h = l), (f == 1 || f == 0) && r != `hsv` && (m = d))
                : ((h = c), (p == 1 || p == 0) && r != `hsv` && (m = u)),
            m === void 0 && (m = u + n * (d - u)),
            (g = f + n * (p - f)),
            r === `oklch` ? new ci([g, m, h], r) : new ci([h, m, g], r)
          );
        },
        ui = li,
        di = function (e, t, n) {
          return ui(e, t, n, `lch`);
        };
      ((Y.lch = di), (Y.hcl = di));
      var fi = b;
      Y.num = function (e, t, n) {
        var r = e.num();
        return new fi(r + n * (t.num() - r), `num`);
      };
      var pi = li;
      Y.hcg = function (e, t, n) {
        return pi(e, t, n, `hcg`);
      };
      var mi = li;
      Y.hsi = function (e, t, n) {
        return mi(e, t, n, `hsi`);
      };
      var hi = li;
      Y.hsl = function (e, t, n) {
        return hi(e, t, n, `hsl`);
      };
      var gi = li;
      Y.hsv = function (e, t, n) {
        return gi(e, t, n, `hsv`);
      };
      var _i = b;
      Y.oklab = function (e, t, n) {
        var r = e.oklab(),
          i = t.oklab();
        return new _i(
          r[0] + n * (i[0] - r[0]),
          r[1] + n * (i[1] - r[1]),
          r[2] + n * (i[2] - r[2]),
          `oklab`,
        );
      };
      var vi = li;
      Y.oklch = function (e, t, n) {
        return vi(e, t, n, `oklch`);
      };
      var yi = b,
        bi = p.clip_rgb,
        xi = Math.pow,
        Si = Math.sqrt,
        Ci = Math.PI,
        wi = Math.cos,
        Ti = Math.sin,
        Ei = Math.atan2,
        Di = function (e, t, n) {
          (t === void 0 && (t = `lrgb`), n === void 0 && (n = null));
          var r = e.length;
          n ||= Array.from(Array(r)).map(function () {
            return 1;
          });
          var i =
            r /
            n.reduce(function (e, t) {
              return e + t;
            });
          if (
            (n.forEach(function (e, t) {
              n[t] *= i;
            }),
            (e = e.map(function (e) {
              return new yi(e);
            })),
            t === `lrgb`)
          )
            return Oi(e, n);
          for (
            var a = e.shift(), o = a.get(t), s = [], c = 0, l = 0, u = 0;
            u < o.length;
            u++
          )
            if (
              ((o[u] = (o[u] || 0) * n[0]),
              s.push(isNaN(o[u]) ? 0 : n[0]),
              t.charAt(u) === `h` && !isNaN(o[u]))
            ) {
              var d = (o[u] / 180) * Ci;
              ((c += wi(d) * n[0]), (l += Ti(d) * n[0]));
            }
          var f = a.alpha() * n[0];
          e.forEach(function (e, r) {
            var i = e.get(t);
            f += e.alpha() * n[r + 1];
            for (var a = 0; a < o.length; a++)
              if (!isNaN(i[a])) {
                if (((s[a] += n[r + 1]), t.charAt(a) === `h`)) {
                  var u = (i[a] / 180) * Ci;
                  ((c += wi(u) * n[r + 1]), (l += Ti(u) * n[r + 1]));
                } else o[a] += i[a] * n[r + 1];
              }
          });
          for (var p = 0; p < o.length; p++)
            if (t.charAt(p) === `h`) {
              for (var m = (Ei(l / s[p], c / s[p]) / Ci) * 180; m < 0;)
                m += 360;
              for (; m >= 360;) m -= 360;
              o[p] = m;
            } else o[p] = o[p] / s[p];
          return ((f /= r), new yi(o, t).alpha(f > 0.99999 ? 1 : f, !0));
        },
        Oi = function (e, t) {
          for (var n = e.length, r = [0, 0, 0, 0], i = 0; i < e.length; i++) {
            var a = e[i],
              o = t[i] / n,
              s = a._rgb;
            ((r[0] += xi(s[0], 2) * o),
              (r[1] += xi(s[1], 2) * o),
              (r[2] += xi(s[2], 2) * o),
              (r[3] += s[3] * o));
          }
          return (
            (r[0] = Si(r[0])),
            (r[1] = Si(r[1])),
            (r[2] = Si(r[2])),
            r[3] > 0.9999999 && (r[3] = 1),
            new yi(bi(r))
          );
        },
        X = S,
        ki = p.type,
        Ai = Math.pow,
        ji = function (e) {
          var t = `rgb`,
            n = X(`#ccc`),
            r = 0,
            i = [0, 1],
            a = [],
            o = [0, 0],
            s = !1,
            c = [],
            l = !1,
            u = 0,
            d = 1,
            f = !1,
            p = {},
            m = !0,
            h = 1,
            g = function (e) {
              if (
                ((e ||= [`#fff`, `#000`]),
                e &&
                  ki(e) === `string` &&
                  X.brewer &&
                  X.brewer[e.toLowerCase()] &&
                  (e = X.brewer[e.toLowerCase()]),
                ki(e) === `array`)
              ) {
                (e.length === 1 && (e = [e[0], e[0]]), (e = e.slice(0)));
                for (var t = 0; t < e.length; t++) e[t] = X(e[t]);
                a.length = 0;
                for (var n = 0; n < e.length; n++) a.push(n / (e.length - 1));
              }
              return (x(), (c = e));
            },
            _ = function (e) {
              if (s != null) {
                for (var t = s.length - 1, n = 0; n < t && e >= s[n];) n++;
                return n - 1;
              }
              return 0;
            },
            v = function (e) {
              return e;
            },
            y = function (e) {
              return e;
            },
            b = function (e, r) {
              var i, l;
              if (((r ??= !1), isNaN(e) || e === null)) return n;
              ((l = r
                ? e
                : s && s.length > 2
                  ? _(e) / (s.length - 2)
                  : d === u
                    ? 1
                    : (e - u) / (d - u)),
                (l = y(l)),
                r || (l = v(l)),
                h !== 1 && (l = Ai(l, h)),
                (l = o[0] + l * (1 - o[0] - o[1])),
                (l = Math.min(1, Math.max(0, l))));
              var f = Math.floor(l * 1e4);
              if (m && p[f]) i = p[f];
              else {
                if (ki(c) === `array`)
                  for (var g = 0; g < a.length; g++) {
                    var b = a[g];
                    if (l <= b) {
                      i = c[g];
                      break;
                    }
                    if (l >= b && g === a.length - 1) {
                      i = c[g];
                      break;
                    }
                    if (l > b && l < a[g + 1]) {
                      ((l = (l - b) / (a[g + 1] - b)),
                        (i = X.interpolate(c[g], c[g + 1], l, t)));
                      break;
                    }
                  }
                else ki(c) === `function` && (i = c(l));
                m && (p[f] = i);
              }
              return i;
            },
            x = function () {
              return (p = {});
            };
          g(e);
          var S = function (e) {
            var t = X(b(e));
            return l && t[l] ? t[l]() : t;
          };
          return (
            (S.classes = function (e) {
              if (e != null) {
                if (ki(e) === `array`) ((s = e), (i = [e[0], e[e.length - 1]]));
                else {
                  var t = X.analyze(i);
                  s = e === 0 ? [t.min, t.max] : X.limits(t, `e`, e);
                }
                return S;
              }
              return s;
            }),
            (S.domain = function (e) {
              if (!arguments.length) return i;
              ((u = e[0]), (d = e[e.length - 1]), (a = []));
              var t = c.length;
              if (e.length === t && u !== d)
                for (var n = 0, r = Array.from(e); n < r.length; n += 1) {
                  var o = r[n];
                  a.push((o - u) / (d - u));
                }
              else {
                for (var s = 0; s < t; s++) a.push(s / (t - 1));
                if (e.length > 2) {
                  var l = e.map(function (t, n) {
                      return n / (e.length - 1);
                    }),
                    f = e.map(function (e) {
                      return (e - u) / (d - u);
                    });
                  f.every(function (e, t) {
                    return l[t] === e;
                  }) ||
                    (y = function (e) {
                      if (e <= 0 || e >= 1) return e;
                      for (var t = 0; e >= f[t + 1];) t++;
                      var n = (e - f[t]) / (f[t + 1] - f[t]);
                      return l[t] + n * (l[t + 1] - l[t]);
                    });
                }
              }
              return ((i = [u, d]), S);
            }),
            (S.mode = function (e) {
              return arguments.length ? ((t = e), x(), S) : t;
            }),
            (S.range = function (e, t) {
              return (g(e), S);
            }),
            (S.out = function (e) {
              return ((l = e), S);
            }),
            (S.spread = function (e) {
              return arguments.length ? ((r = e), S) : r;
            }),
            (S.correctLightness = function (e) {
              return (
                (e ??= !0),
                (f = e),
                x(),
                (v = f
                  ? function (e) {
                      for (
                        var t = b(0, !0).lab()[0],
                          n = b(1, !0).lab()[0],
                          r = t > n,
                          i = b(e, !0).lab()[0],
                          a = t + (n - t) * e,
                          o = i - a,
                          s = 0,
                          c = 1,
                          l = 20;
                        Math.abs(o) > 0.01 && l-- > 0;
                      )
                        (function () {
                          return (
                            r && (o *= -1),
                            o < 0
                              ? ((s = e), (e += (c - e) * 0.5))
                              : ((c = e), (e += (s - e) * 0.5)),
                            (i = b(e, !0).lab()[0]),
                            (o = i - a)
                          );
                        })();
                      return e;
                    }
                  : function (e) {
                      return e;
                    }),
                S
              );
            }),
            (S.padding = function (e) {
              return e == null
                ? o
                : (ki(e) === `number` && (e = [e, e]), (o = e), S);
            }),
            (S.colors = function (t, n) {
              arguments.length < 2 && (n = `hex`);
              var r = [];
              if (arguments.length === 0) r = c.slice(0);
              else if (t === 1) r = [S(0.5)];
              else if (t > 1) {
                var a = i[0],
                  o = i[1] - a;
                r = Mi(0, t, !1).map(function (e) {
                  return S(a + (e / (t - 1)) * o);
                });
              } else {
                e = [];
                var l = [];
                if (s && s.length > 2)
                  for (
                    var u = 1, d = s.length, f = 1 <= d;
                    f ? u < d : u > d;
                    f ? u++ : u--
                  )
                    l.push((s[u - 1] + s[u]) * 0.5);
                else l = i;
                r = l.map(function (e) {
                  return S(e);
                });
              }
              return (
                X[n] &&
                  (r = r.map(function (e) {
                    return e[n]();
                  })),
                r
              );
            }),
            (S.cache = function (e) {
              return e == null ? m : ((m = e), S);
            }),
            (S.gamma = function (e) {
              return e == null ? h : ((h = e), S);
            }),
            (S.nodata = function (e) {
              return e == null ? n : ((n = X(e)), S);
            }),
            S
          );
        };
      function Mi(e, t, n) {
        for (
          var r = [], i = e < t, a = n ? (i ? t + 1 : t - 1) : t, o = e;
          i ? o < a : o > a;
          i ? o++ : o--
        )
          r.push(o);
        return r;
      }
      var Ni = b,
        Pi = ji,
        Fi = function (e) {
          for (var t = [1, 1], n = 1; n < e; n++) {
            for (var r = [1], i = 1; i <= t.length; i++)
              r[i] = (t[i] || 0) + t[i - 1];
            t = r;
          }
          return t;
        },
        Ii = function (e) {
          var t, n, r, i, a, o, s;
          if (
            ((e = e.map(function (e) {
              return new Ni(e);
            })),
            e.length === 2)
          )
            ((t = e.map(function (e) {
              return e.lab();
            })),
              (a = t[0]),
              (o = t[1]),
              (i = function (e) {
                return new Ni(
                  [0, 1, 2].map(function (t) {
                    return a[t] + e * (o[t] - a[t]);
                  }),
                  `lab`,
                );
              }));
          else if (e.length === 3)
            ((n = e.map(function (e) {
              return e.lab();
            })),
              (a = n[0]),
              (o = n[1]),
              (s = n[2]),
              (i = function (e) {
                return new Ni(
                  [0, 1, 2].map(function (t) {
                    return (
                      (1 - e) * (1 - e) * a[t] +
                      2 * (1 - e) * e * o[t] +
                      e * e * s[t]
                    );
                  }),
                  `lab`,
                );
              }));
          else if (e.length === 4) {
            var c;
            ((r = e.map(function (e) {
              return e.lab();
            })),
              (a = r[0]),
              (o = r[1]),
              (s = r[2]),
              (c = r[3]),
              (i = function (e) {
                return new Ni(
                  [0, 1, 2].map(function (t) {
                    return (
                      (1 - e) * (1 - e) * (1 - e) * a[t] +
                      3 * (1 - e) * (1 - e) * e * o[t] +
                      3 * (1 - e) * e * e * s[t] +
                      e * e * e * c[t]
                    );
                  }),
                  `lab`,
                );
              }));
          } else if (e.length >= 5) {
            var l = e.map(function (e) {
                return e.lab();
              }),
              u,
              d = e.length - 1;
            ((u = Fi(d)),
              (i = function (e) {
                var t = 1 - e;
                return new Ni(
                  [0, 1, 2].map(function (n) {
                    return l.reduce(function (r, i, a) {
                      return r + u[a] * t ** +(d - a) * e ** +a * i[n];
                    }, 0);
                  }),
                  `lab`,
                );
              }));
          } else
            throw RangeError(`No point in running bezier with only one color.`);
          return i;
        },
        Li = function (e) {
          var t = Ii(e);
          return (
            (t.scale = function () {
              return Pi(t);
            }),
            t
          );
        },
        Ri = S,
        Z = function (e, t, n) {
          if (!Z[n]) throw Error(`unknown blend mode ` + n);
          return Z[n](e, t);
        },
        zi = function (e) {
          return function (t, n) {
            var r = Ri(n).rgb(),
              i = Ri(t).rgb();
            return Ri.rgb(e(r, i));
          };
        },
        Bi = function (e) {
          return function (t, n) {
            var r = [];
            return (
              (r[0] = e(t[0], n[0])),
              (r[1] = e(t[1], n[1])),
              (r[2] = e(t[2], n[2])),
              r
            );
          };
        };
      ((Z.normal = zi(
        Bi(function (e) {
          return e;
        }),
      )),
        (Z.multiply = zi(
          Bi(function (e, t) {
            return (e * t) / 255;
          }),
        )),
        (Z.screen = zi(
          Bi(function (e, t) {
            return 255 * (1 - (1 - e / 255) * (1 - t / 255));
          }),
        )),
        (Z.overlay = zi(
          Bi(function (e, t) {
            return t < 128
              ? (2 * e * t) / 255
              : 255 * (1 - 2 * (1 - e / 255) * (1 - t / 255));
          }),
        )),
        (Z.darken = zi(
          Bi(function (e, t) {
            return e > t ? t : e;
          }),
        )),
        (Z.lighten = zi(
          Bi(function (e, t) {
            return e > t ? e : t;
          }),
        )),
        (Z.dodge = zi(
          Bi(function (e, t) {
            return e === 255
              ? 255
              : ((e = ((t / 255) * 255) / (1 - e / 255)), e > 255 ? 255 : e);
          }),
        )),
        (Z.burn = zi(
          Bi(function (e, t) {
            return 255 * (1 - (1 - t / 255) / (e / 255));
          }),
        )));
      for (
        var Vi = Z,
          Hi = p.type,
          Ui = p.clip_rgb,
          Wi = p.TWOPI,
          Gi = Math.pow,
          Ki = Math.sin,
          qi = Math.cos,
          Ji = S,
          Yi = function (e, t, n, r, i) {
            (e === void 0 && (e = 300),
              t === void 0 && (t = -1.5),
              n === void 0 && (n = 1),
              r === void 0 && (r = 1),
              i === void 0 && (i = [0, 1]));
            var a = 0,
              o;
            Hi(i) === `array` ? (o = i[1] - i[0]) : ((o = 0), (i = [i, i]));
            var s = function (s) {
              var c = Wi * ((e + 120) / 360 + t * s),
                l = Gi(i[0] + o * s, r),
                u = ((a === 0 ? n : n[0] + s * a) * l * (1 - l)) / 2,
                d = qi(c),
                f = Ki(c),
                p = l + u * (-0.14861 * d + 1.78277 * f),
                m = l + u * (-0.29227 * d - 0.90649 * f),
                h = l + 1.97294 * d * u;
              return Ji(Ui([p * 255, m * 255, h * 255, 1]));
            };
            return (
              (s.start = function (t) {
                return t == null ? e : ((e = t), s);
              }),
              (s.rotations = function (e) {
                return e == null ? t : ((t = e), s);
              }),
              (s.gamma = function (e) {
                return e == null ? r : ((r = e), s);
              }),
              (s.hue = function (e) {
                return e == null
                  ? n
                  : ((n = e),
                    Hi(n) === `array`
                      ? ((a = n[1] - n[0]), a === 0 && (n = n[1]))
                      : (a = 0),
                    s);
              }),
              (s.lightness = function (e) {
                return e == null
                  ? i
                  : (Hi(e) === `array`
                      ? ((i = e), (o = e[1] - e[0]))
                      : ((i = [e, e]), (o = 0)),
                    s);
              }),
              (s.scale = function () {
                return Ji.scale(s);
              }),
              s.hue(n),
              s
            );
          },
          Xi = b,
          Zi = `0123456789abcdef`,
          Qi = Math.floor,
          $i = Math.random,
          ea = function () {
            for (var e = `#`, t = 0; t < 6; t++) e += Zi.charAt(Qi($i() * 16));
            return new Xi(e, `hex`);
          },
          ta = s,
          na = Math.log,
          ra = Math.pow,
          ia = Math.floor,
          aa = Math.abs,
          oa = function (e, t) {
            t === void 0 && (t = null);
            var n = {
              min: Number.MAX_VALUE,
              max: Number.MAX_VALUE * -1,
              sum: 0,
              values: [],
              count: 0,
            };
            return (
              ta(e) === `object` && (e = Object.values(e)),
              e.forEach(function (e) {
                (t && ta(e) === `object` && (e = e[t]),
                  e != null &&
                    !isNaN(e) &&
                    (n.values.push(e),
                    (n.sum += e),
                    e < n.min && (n.min = e),
                    e > n.max && (n.max = e),
                    (n.count += 1)));
              }),
              (n.domain = [n.min, n.max]),
              (n.limits = function (e, t) {
                return sa(n, e, t);
              }),
              n
            );
          },
          sa = function (e, t, n) {
            (t === void 0 && (t = `equal`),
              n === void 0 && (n = 7),
              ta(e) == `array` && (e = oa(e)));
            var r = e.min,
              i = e.max,
              a = e.values.sort(function (e, t) {
                return e - t;
              });
            if (n === 1) return [r, i];
            var o = [];
            if (
              (t.substr(0, 1) === `c` && (o.push(r), o.push(i)),
              t.substr(0, 1) === `e`)
            ) {
              o.push(r);
              for (var s = 1; s < n; s++) o.push(r + (s / n) * (i - r));
              o.push(i);
            } else if (t.substr(0, 1) === `l`) {
              if (r <= 0)
                throw Error(
                  `Logarithmic scales are only possible for values > 0`,
                );
              var c = Math.LOG10E * na(r),
                l = Math.LOG10E * na(i);
              o.push(r);
              for (var u = 1; u < n; u++) o.push(ra(10, c + (u / n) * (l - c)));
              o.push(i);
            } else if (t.substr(0, 1) === `q`) {
              o.push(r);
              for (var d = 1; d < n; d++) {
                var f = ((a.length - 1) * d) / n,
                  p = ia(f);
                if (p === f) o.push(a[p]);
                else {
                  var m = f - p;
                  o.push(a[p] * (1 - m) + a[p + 1] * m);
                }
              }
              o.push(i);
            } else if (t.substr(0, 1) === `k`) {
              var h,
                g = a.length,
                _ = Array(g),
                v = Array(n),
                y = !0,
                b = 0,
                x = null;
              ((x = []), x.push(r));
              for (var S = 1; S < n; S++) x.push(r + (S / n) * (i - r));
              for (x.push(i); y;) {
                for (var C = 0; C < n; C++) v[C] = 0;
                for (var w = 0; w < g; w++)
                  for (
                    var T = a[w], E = Number.MAX_VALUE, D = void 0, O = 0;
                    O < n;
                    O++
                  ) {
                    var k = aa(x[O] - T);
                    (k < E && ((E = k), (D = O)), v[D]++, (_[w] = D));
                  }
                for (var A = Array(n), ee = 0; ee < n; ee++) A[ee] = null;
                for (var j = 0; j < g; j++)
                  ((h = _[j]), A[h] === null ? (A[h] = a[j]) : (A[h] += a[j]));
                for (var te = 0; te < n; te++) A[te] *= 1 / v[te];
                y = !1;
                for (var M = 0; M < n; M++)
                  if (A[M] !== x[M]) {
                    y = !0;
                    break;
                  }
                ((x = A), b++, b > 200 && (y = !1));
              }
              for (var N = {}, P = 0; P < n; P++) N[P] = [];
              for (var ne = 0; ne < g; ne++) ((h = _[ne]), N[h].push(a[ne]));
              for (var F = [], I = 0; I < n; I++)
                (F.push(N[I][0]), F.push(N[I][N[I].length - 1]));
              ((F = F.sort(function (e, t) {
                return e - t;
              })),
                o.push(F[0]));
              for (var re = 1; re < F.length; re += 2) {
                var ie = F[re];
                !isNaN(ie) && o.indexOf(ie) === -1 && o.push(ie);
              }
            }
            return o;
          },
          ca = { analyze: oa, limits: sa },
          la = b,
          ua = function (e, t) {
            ((e = new la(e)), (t = new la(t)));
            var n = e.luminance(),
              r = t.luminance();
            return n > r ? (n + 0.05) / (r + 0.05) : (r + 0.05) / (n + 0.05);
          },
          da = b,
          fa = Math.sqrt,
          Q = Math.pow,
          pa = Math.min,
          ma = Math.max,
          ha = Math.atan2,
          ga = Math.abs,
          _a = Math.cos,
          va = Math.sin,
          ya = Math.exp,
          ba = Math.PI,
          xa = function (e, t, n, r, i) {
            (n === void 0 && (n = 1),
              r === void 0 && (r = 1),
              i === void 0 && (i = 1));
            var a = function (e) {
                return (360 * e) / (2 * ba);
              },
              o = function (e) {
                return (2 * ba * e) / 360;
              };
            ((e = new da(e)), (t = new da(t)));
            var s = Array.from(e.lab()),
              c = s[0],
              l = s[1],
              u = s[2],
              d = Array.from(t.lab()),
              f = d[0],
              p = d[1],
              m = d[2],
              h = (c + f) / 2,
              g = (fa(Q(l, 2) + Q(u, 2)) + fa(Q(p, 2) + Q(m, 2))) / 2,
              _ = 0.5 * (1 - fa(Q(g, 7) / (Q(g, 7) + Q(25, 7)))),
              v = l * (1 + _),
              y = p * (1 + _),
              b = fa(Q(v, 2) + Q(u, 2)),
              x = fa(Q(y, 2) + Q(m, 2)),
              S = (b + x) / 2,
              C = a(ha(u, v)),
              w = a(ha(m, y)),
              T = C >= 0 ? C : C + 360,
              E = w >= 0 ? w : w + 360,
              D = ga(T - E) > 180 ? (T + E + 360) / 2 : (T + E) / 2,
              O =
                1 -
                0.17 * _a(o(D - 30)) +
                0.24 * _a(o(2 * D)) +
                0.32 * _a(o(3 * D + 6)) -
                0.2 * _a(o(4 * D - 63)),
              k = E - T;
            ((k = ga(k) <= 180 ? k : E <= T ? k + 360 : k - 360),
              (k = 2 * fa(b * x) * va(o(k) / 2)));
            var A = f - c,
              ee = x - b,
              j = 1 + (0.015 * Q(h - 50, 2)) / fa(20 + Q(h - 50, 2)),
              te = 1 + 0.045 * S,
              M = 1 + 0.015 * S * O,
              N = 30 * ya(-Q((D - 275) / 25, 2)),
              P = -(2 * fa(Q(S, 7) / (Q(S, 7) + Q(25, 7)))) * va(2 * o(N));
            return ma(
              0,
              pa(
                100,
                fa(
                  Q(A / (n * j), 2) +
                    Q(ee / (r * te), 2) +
                    Q(k / (i * M), 2) +
                    P * (ee / (r * te)) * (k / (i * M)),
                ),
              ),
            );
          },
          Sa = b,
          Ca = function (e, t, n) {
            (n === void 0 && (n = `lab`), (e = new Sa(e)), (t = new Sa(t)));
            var r = e.get(n),
              i = t.get(n),
              a = 0;
            for (var o in r) {
              var s = (r[o] || 0) - (i[o] || 0);
              a += s * s;
            }
            return Math.sqrt(a);
          },
          wa = b,
          Ta = function () {
            for (var e = [], t = arguments.length; t--;) e[t] = arguments[t];
            try {
              return (
                new (Function.prototype.bind.apply(wa, [null].concat(e)))(),
                !0
              );
            } catch {
              return !1;
            }
          },
          Ea = S,
          Da = ji,
          Oa = {
            cool: function () {
              return Da([Ea.hsl(180, 1, 0.9), Ea.hsl(250, 0.7, 0.4)]);
            },
            hot: function () {
              return Da([`#000`, `#f00`, `#ff0`, `#fff`]).mode(`rgb`);
            },
          },
          ka = {
            OrRd: [
              `#fff7ec`,
              `#fee8c8`,
              `#fdd49e`,
              `#fdbb84`,
              `#fc8d59`,
              `#ef6548`,
              `#d7301f`,
              `#b30000`,
              `#7f0000`,
            ],
            PuBu: [
              `#fff7fb`,
              `#ece7f2`,
              `#d0d1e6`,
              `#a6bddb`,
              `#74a9cf`,
              `#3690c0`,
              `#0570b0`,
              `#045a8d`,
              `#023858`,
            ],
            BuPu: [
              `#f7fcfd`,
              `#e0ecf4`,
              `#bfd3e6`,
              `#9ebcda`,
              `#8c96c6`,
              `#8c6bb1`,
              `#88419d`,
              `#810f7c`,
              `#4d004b`,
            ],
            Oranges: [
              `#fff5eb`,
              `#fee6ce`,
              `#fdd0a2`,
              `#fdae6b`,
              `#fd8d3c`,
              `#f16913`,
              `#d94801`,
              `#a63603`,
              `#7f2704`,
            ],
            BuGn: [
              `#f7fcfd`,
              `#e5f5f9`,
              `#ccece6`,
              `#99d8c9`,
              `#66c2a4`,
              `#41ae76`,
              `#238b45`,
              `#006d2c`,
              `#00441b`,
            ],
            YlOrBr: [
              `#ffffe5`,
              `#fff7bc`,
              `#fee391`,
              `#fec44f`,
              `#fe9929`,
              `#ec7014`,
              `#cc4c02`,
              `#993404`,
              `#662506`,
            ],
            YlGn: [
              `#ffffe5`,
              `#f7fcb9`,
              `#d9f0a3`,
              `#addd8e`,
              `#78c679`,
              `#41ab5d`,
              `#238443`,
              `#006837`,
              `#004529`,
            ],
            Reds: [
              `#fff5f0`,
              `#fee0d2`,
              `#fcbba1`,
              `#fc9272`,
              `#fb6a4a`,
              `#ef3b2c`,
              `#cb181d`,
              `#a50f15`,
              `#67000d`,
            ],
            RdPu: [
              `#fff7f3`,
              `#fde0dd`,
              `#fcc5c0`,
              `#fa9fb5`,
              `#f768a1`,
              `#dd3497`,
              `#ae017e`,
              `#7a0177`,
              `#49006a`,
            ],
            Greens: [
              `#f7fcf5`,
              `#e5f5e0`,
              `#c7e9c0`,
              `#a1d99b`,
              `#74c476`,
              `#41ab5d`,
              `#238b45`,
              `#006d2c`,
              `#00441b`,
            ],
            YlGnBu: [
              `#ffffd9`,
              `#edf8b1`,
              `#c7e9b4`,
              `#7fcdbb`,
              `#41b6c4`,
              `#1d91c0`,
              `#225ea8`,
              `#253494`,
              `#081d58`,
            ],
            Purples: [
              `#fcfbfd`,
              `#efedf5`,
              `#dadaeb`,
              `#bcbddc`,
              `#9e9ac8`,
              `#807dba`,
              `#6a51a3`,
              `#54278f`,
              `#3f007d`,
            ],
            GnBu: [
              `#f7fcf0`,
              `#e0f3db`,
              `#ccebc5`,
              `#a8ddb5`,
              `#7bccc4`,
              `#4eb3d3`,
              `#2b8cbe`,
              `#0868ac`,
              `#084081`,
            ],
            Greys: [
              `#ffffff`,
              `#f0f0f0`,
              `#d9d9d9`,
              `#bdbdbd`,
              `#969696`,
              `#737373`,
              `#525252`,
              `#252525`,
              `#000000`,
            ],
            YlOrRd: [
              `#ffffcc`,
              `#ffeda0`,
              `#fed976`,
              `#feb24c`,
              `#fd8d3c`,
              `#fc4e2a`,
              `#e31a1c`,
              `#bd0026`,
              `#800026`,
            ],
            PuRd: [
              `#f7f4f9`,
              `#e7e1ef`,
              `#d4b9da`,
              `#c994c7`,
              `#df65b0`,
              `#e7298a`,
              `#ce1256`,
              `#980043`,
              `#67001f`,
            ],
            Blues: [
              `#f7fbff`,
              `#deebf7`,
              `#c6dbef`,
              `#9ecae1`,
              `#6baed6`,
              `#4292c6`,
              `#2171b5`,
              `#08519c`,
              `#08306b`,
            ],
            PuBuGn: [
              `#fff7fb`,
              `#ece2f0`,
              `#d0d1e6`,
              `#a6bddb`,
              `#67a9cf`,
              `#3690c0`,
              `#02818a`,
              `#016c59`,
              `#014636`,
            ],
            Viridis: [
              `#440154`,
              `#482777`,
              `#3f4a8a`,
              `#31678e`,
              `#26838f`,
              `#1f9d8a`,
              `#6cce5a`,
              `#b6de2b`,
              `#fee825`,
            ],
            Spectral: [
              `#9e0142`,
              `#d53e4f`,
              `#f46d43`,
              `#fdae61`,
              `#fee08b`,
              `#ffffbf`,
              `#e6f598`,
              `#abdda4`,
              `#66c2a5`,
              `#3288bd`,
              `#5e4fa2`,
            ],
            RdYlGn: [
              `#a50026`,
              `#d73027`,
              `#f46d43`,
              `#fdae61`,
              `#fee08b`,
              `#ffffbf`,
              `#d9ef8b`,
              `#a6d96a`,
              `#66bd63`,
              `#1a9850`,
              `#006837`,
            ],
            RdBu: [
              `#67001f`,
              `#b2182b`,
              `#d6604d`,
              `#f4a582`,
              `#fddbc7`,
              `#f7f7f7`,
              `#d1e5f0`,
              `#92c5de`,
              `#4393c3`,
              `#2166ac`,
              `#053061`,
            ],
            PiYG: [
              `#8e0152`,
              `#c51b7d`,
              `#de77ae`,
              `#f1b6da`,
              `#fde0ef`,
              `#f7f7f7`,
              `#e6f5d0`,
              `#b8e186`,
              `#7fbc41`,
              `#4d9221`,
              `#276419`,
            ],
            PRGn: [
              `#40004b`,
              `#762a83`,
              `#9970ab`,
              `#c2a5cf`,
              `#e7d4e8`,
              `#f7f7f7`,
              `#d9f0d3`,
              `#a6dba0`,
              `#5aae61`,
              `#1b7837`,
              `#00441b`,
            ],
            RdYlBu: [
              `#a50026`,
              `#d73027`,
              `#f46d43`,
              `#fdae61`,
              `#fee090`,
              `#ffffbf`,
              `#e0f3f8`,
              `#abd9e9`,
              `#74add1`,
              `#4575b4`,
              `#313695`,
            ],
            BrBG: [
              `#543005`,
              `#8c510a`,
              `#bf812d`,
              `#dfc27d`,
              `#f6e8c3`,
              `#f5f5f5`,
              `#c7eae5`,
              `#80cdc1`,
              `#35978f`,
              `#01665e`,
              `#003c30`,
            ],
            RdGy: [
              `#67001f`,
              `#b2182b`,
              `#d6604d`,
              `#f4a582`,
              `#fddbc7`,
              `#ffffff`,
              `#e0e0e0`,
              `#bababa`,
              `#878787`,
              `#4d4d4d`,
              `#1a1a1a`,
            ],
            PuOr: [
              `#7f3b08`,
              `#b35806`,
              `#e08214`,
              `#fdb863`,
              `#fee0b6`,
              `#f7f7f7`,
              `#d8daeb`,
              `#b2abd2`,
              `#8073ac`,
              `#542788`,
              `#2d004b`,
            ],
            Set2: [
              `#66c2a5`,
              `#fc8d62`,
              `#8da0cb`,
              `#e78ac3`,
              `#a6d854`,
              `#ffd92f`,
              `#e5c494`,
              `#b3b3b3`,
            ],
            Accent: [
              `#7fc97f`,
              `#beaed4`,
              `#fdc086`,
              `#ffff99`,
              `#386cb0`,
              `#f0027f`,
              `#bf5b17`,
              `#666666`,
            ],
            Set1: [
              `#e41a1c`,
              `#377eb8`,
              `#4daf4a`,
              `#984ea3`,
              `#ff7f00`,
              `#ffff33`,
              `#a65628`,
              `#f781bf`,
              `#999999`,
            ],
            Set3: [
              `#8dd3c7`,
              `#ffffb3`,
              `#bebada`,
              `#fb8072`,
              `#80b1d3`,
              `#fdb462`,
              `#b3de69`,
              `#fccde5`,
              `#d9d9d9`,
              `#bc80bd`,
              `#ccebc5`,
              `#ffed6f`,
            ],
            Dark2: [
              `#1b9e77`,
              `#d95f02`,
              `#7570b3`,
              `#e7298a`,
              `#66a61e`,
              `#e6ab02`,
              `#a6761d`,
              `#666666`,
            ],
            Paired: [
              `#a6cee3`,
              `#1f78b4`,
              `#b2df8a`,
              `#33a02c`,
              `#fb9a99`,
              `#e31a1c`,
              `#fdbf6f`,
              `#ff7f00`,
              `#cab2d6`,
              `#6a3d9a`,
              `#ffff99`,
              `#b15928`,
            ],
            Pastel2: [
              `#b3e2cd`,
              `#fdcdac`,
              `#cbd5e8`,
              `#f4cae4`,
              `#e6f5c9`,
              `#fff2ae`,
              `#f1e2cc`,
              `#cccccc`,
            ],
            Pastel1: [
              `#fbb4ae`,
              `#b3cde3`,
              `#ccebc5`,
              `#decbe4`,
              `#fed9a6`,
              `#ffffcc`,
              `#e5d8bd`,
              `#fddaec`,
              `#f2f2f2`,
            ],
          },
          Aa = 0,
          ja = Object.keys(ka);
        Aa < ja.length;
        Aa += 1
      ) {
        var Ma = ja[Aa];
        ka[Ma.toLowerCase()] = ka[Ma];
      }
      var Na = ka,
        $ = S;
      return (
        ($.average = Di),
        ($.bezier = Li),
        ($.blend = Vi),
        ($.cubehelix = Yi),
        ($.mix = $.interpolate = Yr),
        ($.random = ea),
        ($.scale = ji),
        ($.analyze = ca.analyze),
        ($.contrast = ua),
        ($.deltaE = xa),
        ($.distance = Ca),
        ($.limits = ca.limits),
        ($.valid = Ta),
        ($.scales = Oa),
        ($.colors = wn),
        ($.brewer = Na),
        $
      );
    });
  }),
  f = class e {
    static clone(t) {
      let n = Array.isArray(t) ? [] : {};
      for (let r in t) {
        let i = t[r];
        n[r] =
          typeof i.clone == `function`
            ? i.clone()
            : typeof i == `object`
              ? e.clone(i)
              : i;
      }
      return n;
    }
    static equals(e, t) {
      if (e.length !== t.length) return !1;
      let n = e.slice().sort(),
        r = t.slice().sort();
      for (let e = 0; e < n.length; e++) if (n[e] !== r[e]) return !1;
      return !0;
    }
    static print(e) {
      if (e.length == 0) return ``;
      let t = `(`;
      for (let n = 0; n < e.length; n++)
        t += e[n].id ? e[n].id + `, ` : e[n] + `, `;
      return ((t = t.substring(0, t.length - 2)), t + `)`);
    }
    static each(e, t) {
      for (let n = 0; n < e.length; n++) t(e[n]);
    }
    static get(e, t, n) {
      for (let r = 0; r < e.length; r++) if (e[r][t] == n) return e[r];
    }
    static contains(e, t) {
      if (!t.property && !t.func) {
        for (let n = 0; n < e.length; n++) if (e[n] == t.value) return !0;
      } else if (t.func) {
        for (let n = 0; n < e.length; n++) if (t.func(e[n])) return !0;
      } else
        for (let n = 0; n < e.length; n++)
          if (e[n][t.property] == t.value) return !0;
      return !1;
    }
    static intersection(e, t) {
      let n = [];
      for (let r = 0; r < e.length; r++)
        for (let i = 0; i < t.length; i++) e[r] === t[i] && n.push(e[r]);
      return n;
    }
    static unique(e) {
      let t = {};
      return e.filter(function (e) {
        return t[e] === void 0 ? (t[e] = !0) : !1;
      });
    }
    static count(e, t) {
      let n = 0;
      for (let r = 0; r < e.length; r++) e[r] === t && n++;
      return n;
    }
    static toggle(e, t) {
      let n = [],
        r = !1;
      for (let i = 0; i < e.length; i++) e[i] === t ? (r = !0) : n.push(e[i]);
      return (r || n.push(t), n);
    }
    static remove(e, t) {
      let n = [];
      for (let r = 0; r < e.length; r++) e[r] !== t && n.push(e[r]);
      return n;
    }
    static removeUnique(e, t) {
      let n = e.indexOf(t);
      return (n > -1 && e.splice(n, 1), e);
    }
    static removeAll(e, t) {
      return e.filter(function (e) {
        return t.indexOf(e) === -1;
      });
    }
    static merge(e, t) {
      let n = Array(e.length + t.length);
      for (let t = 0; t < e.length; t++) n[t] = e[t];
      for (let r = 0; r < t.length; r++) n[e.length + r] = t[r];
      return n;
    }
    static containsAll(e, t) {
      let n = 0;
      for (let r = 0; r < e.length; r++)
        for (let i = 0; i < t.length; i++) e[r] === t[i] && n++;
      return n === t.length;
    }
    static sortByAtomicNumberDesc(e) {
      let t = e.map(function (e, t) {
        return { index: t, value: e.atomicNumber.split(`.`).map(Number) };
      });
      return (
        t.sort(function (e, t) {
          let n = Math.min(t.value.length, e.value.length),
            r = 0;
          for (; r < n && t.value[r] === e.value[r];) r++;
          return r === n
            ? t.value.length - e.value.length
            : t.value[r] - e.value[r];
        }),
        t.map(function (t) {
          return e[t.index];
        })
      );
    }
    static deepCopy(t) {
      let n = [];
      for (let r = 0; r < t.length; r++) {
        let i = t[r];
        i instanceof Array ? (n[r] = e.deepCopy(i)) : (n[r] = i);
      }
      return n;
    }
  },
  p = class e {
    constructor(e, t = `-`) {
      ((this.idx = null),
        (this.element = e.length === 1 ? e.toUpperCase() : e),
        (this.drawExplicit = !1),
        (this.ringbonds = []),
        (this.rings = []),
        (this.bondType = t),
        (this.branchBond = null),
        (this.isBridge = !1),
        (this.isBridgeNode = !1),
        (this.originalRings = []),
        (this.bridgedRing = null),
        (this.anchoredRings = []),
        (this.bracket = null),
        (this.plane = 0),
        (this.attachedPseudoElements = {}),
        (this.hasAttachedPseudoElements = !1),
        (this.isDrawn = !0),
        (this.isConnectedToRing = !1),
        (this.neighbouringElements = []),
        (this.isPartOfAromaticRing = e !== this.element),
        (this.bondCount = 0),
        (this.chirality = ``),
        (this.isStereoCenter = !1),
        (this.priority = 0),
        (this.mainChain = !1),
        (this.subtreeDepth = 1),
        (this.class = void 0));
    }
    addNeighbouringElement(e) {
      this.neighbouringElements.push(e);
    }
    attachPseudoElement(e, t, n = 0, r = 0) {
      (n === null && (n = 0), r === null && (r = 0));
      let i = n + e + r;
      (this.attachedPseudoElements[i]
        ? (this.attachedPseudoElements[i].count += 1)
        : (this.attachedPseudoElements[i] = {
            element: e,
            count: 1,
            hydrogenCount: n,
            previousElement: t,
            charge: r,
          }),
        (this.hasAttachedPseudoElements = !0));
    }
    getAttachedPseudoElements() {
      let e = {};
      return (
        Object.keys(this.attachedPseudoElements)
          .sort()
          .forEach((t) => {
            e[t] = this.attachedPseudoElements[t];
          }),
        e
      );
    }
    getAttachedPseudoElementsCount() {
      return Object.keys(this.attachedPseudoElements).length;
    }
    isHeteroAtom() {
      return this.element !== `C` && this.element !== `H`;
    }
    addAnchoredRing(e) {
      f.contains(this.anchoredRings, { value: e }) ||
        this.anchoredRings.push(e);
    }
    getRingbondCount() {
      return this.ringbonds.length;
    }
    backupRings() {
      this.originalRings = Array(this.rings.length);
      for (let e = 0; e < this.rings.length; e++)
        this.originalRings[e] = this.rings[e];
    }
    restoreRings() {
      this.rings = Array(this.originalRings.length);
      for (let e = 0; e < this.originalRings.length; e++)
        this.rings[e] = this.originalRings[e];
    }
    haveCommonRingbond(e, t) {
      for (let n = 0; n < e.ringbonds.length; n++)
        for (let r = 0; r < t.ringbonds.length; r++)
          if (e.ringbonds[n].id == t.ringbonds[r].id) return !0;
      return !1;
    }
    neighbouringElementsEqual(e) {
      if (e.length !== this.neighbouringElements.length) return !1;
      (e.sort(), this.neighbouringElements.sort());
      for (let t = 0; t < this.neighbouringElements.length; t++)
        if (e[t] !== this.neighbouringElements[t]) return !1;
      return !0;
    }
    getAtomicNumber() {
      return e.atomicNumbers[this.element];
    }
    countImplicitHydrogens() {
      if (this.bracket)
        return this.bracket.chirality ? 0 : this.bracket.hcount || 0;
      let t = this.bondCount;
      if (this.isPartOfAromaticRing) {
        if (this.element !== `C`) return 0;
        t += 1;
      }
      let n = e.VALENCES[this.element];
      if (n === void 0) return 0;
      let r = n.find((e) => e >= t);
      return r === void 0 ? 0 : r - t;
    }
    static get atomicNumbers() {
      return {
        H: 1,
        He: 2,
        Li: 3,
        Be: 4,
        B: 5,
        b: 5,
        C: 6,
        c: 6,
        N: 7,
        n: 7,
        O: 8,
        o: 8,
        F: 9,
        Ne: 10,
        Na: 11,
        Mg: 12,
        Al: 13,
        Si: 14,
        P: 15,
        p: 15,
        S: 16,
        s: 16,
        Cl: 17,
        Ar: 18,
        K: 19,
        Ca: 20,
        Sc: 21,
        Ti: 22,
        V: 23,
        Cr: 24,
        Mn: 25,
        Fe: 26,
        Co: 27,
        Ni: 28,
        Cu: 29,
        Zn: 30,
        Ga: 31,
        Ge: 32,
        As: 33,
        Se: 34,
        Br: 35,
        Kr: 36,
        Rb: 37,
        Sr: 38,
        Y: 39,
        Zr: 40,
        Nb: 41,
        Mo: 42,
        Tc: 43,
        Ru: 44,
        Rh: 45,
        Pd: 46,
        Ag: 47,
        Cd: 48,
        In: 49,
        Sn: 50,
        Sb: 51,
        Te: 52,
        I: 53,
        Xe: 54,
        Cs: 55,
        Ba: 56,
        La: 57,
        Ce: 58,
        Pr: 59,
        Nd: 60,
        Pm: 61,
        Sm: 62,
        Eu: 63,
        Gd: 64,
        Tb: 65,
        Dy: 66,
        Ho: 67,
        Er: 68,
        Tm: 69,
        Yb: 70,
        Lu: 71,
        Hf: 72,
        Ta: 73,
        W: 74,
        Re: 75,
        Os: 76,
        Ir: 77,
        Pt: 78,
        Au: 79,
        Hg: 80,
        Tl: 81,
        Pb: 82,
        Bi: 83,
        Po: 84,
        At: 85,
        Rn: 86,
        Fr: 87,
        Ra: 88,
        Ac: 89,
        Th: 90,
        Pa: 91,
        U: 92,
        Np: 93,
        Pu: 94,
        Am: 95,
        Cm: 96,
        Bk: 97,
        Cf: 98,
        Es: 99,
        Fm: 100,
        Md: 101,
        No: 102,
        Lr: 103,
        Rf: 104,
        Db: 105,
        Sg: 106,
        Bh: 107,
        Hs: 108,
        Mt: 109,
        Ds: 110,
        Rg: 111,
        Cn: 112,
        Uut: 113,
        Uuq: 114,
        Uup: 115,
        Uuh: 116,
        Uus: 117,
        Uuo: 118,
      };
    }
  };
u(p, `VALENCES`, {
  H: [1],
  B: [3],
  C: [4],
  N: [3, 5],
  O: [2],
  F: [1],
  P: [3, 5],
  S: [2, 4, 6],
  Cl: [1],
  Br: [1],
  I: [1],
});
var m = p,
  h = class e {
    constructor(t, n) {
      arguments.length == 0
        ? ((this.x = 0), (this.y = 0))
        : t instanceof e
          ? ((this.x = t.x), (this.y = t.y))
          : ((this.x = t), (this.y = n));
    }
    clone() {
      return new e(this.x, this.y);
    }
    toString() {
      return `(` + this.x + `,` + this.y + `)`;
    }
    add(e) {
      return ((this.x += e.x), (this.y += e.y), this);
    }
    subtract(e) {
      return ((this.x -= e.x), (this.y -= e.y), this);
    }
    divide(e) {
      return ((this.x /= e), (this.y /= e), this);
    }
    multiply(e) {
      return ((this.x *= e.x), (this.y *= e.y), this);
    }
    multiplyScalar(e) {
      return ((this.x *= e), (this.y *= e), this);
    }
    invert() {
      return ((this.x = -this.x), (this.y = -this.y), this);
    }
    angle() {
      return Math.atan2(this.y, this.x);
    }
    distance(e) {
      return Math.sqrt(
        (e.x - this.x) * (e.x - this.x) + (e.y - this.y) * (e.y - this.y),
      );
    }
    distanceSq(e) {
      return (e.x - this.x) * (e.x - this.x) + (e.y - this.y) * (e.y - this.y);
    }
    clockwise(e) {
      let t = this.y * e.x,
        n = this.x * e.y;
      return t > n ? -1 : t === n ? 0 : 1;
    }
    relativeClockwise(e, t) {
      let n = (this.y - e.y) * (t.x - e.x),
        r = (this.x - e.x) * (t.y - e.y);
      return n > r ? -1 : n === r ? 0 : 1;
    }
    rotate(t) {
      let n = new e(0, 0),
        r = Math.cos(t),
        i = Math.sin(t);
      return (
        (n.x = this.x * r - this.y * i),
        (n.y = this.x * i + this.y * r),
        (this.x = n.x),
        (this.y = n.y),
        this
      );
    }
    rotateAround(e, t) {
      let n = Math.sin(e),
        r = Math.cos(e);
      ((this.x -= t.x), (this.y -= t.y));
      let i = this.x * r - this.y * n,
        a = this.x * n + this.y * r;
      return ((this.x = i + t.x), (this.y = a + t.y), this);
    }
    rotateTo(t, n, r = 0) {
      ((this.x += 0.001), (this.y -= 0.001));
      let i = e.subtract(this, n),
        a = e.subtract(t, n),
        o = e.angle(a, i);
      return (this.rotateAround(o + r, n), this);
    }
    rotateAwayFrom(e, t, n) {
      this.rotateAround(n, t);
      let r = this.distanceSq(e);
      (this.rotateAround(-2 * n, t),
        this.distanceSq(e) < r && this.rotateAround(2 * n, t));
    }
    getRotateAwayFromAngle(e, t, n) {
      let r = this.clone();
      r.rotateAround(n, t);
      let i = r.distanceSq(e);
      return (r.rotateAround(-2 * n, t), r.distanceSq(e) < i ? n : -n);
    }
    getRotateTowardsAngle(e, t, n) {
      let r = this.clone();
      r.rotateAround(n, t);
      let i = r.distanceSq(e);
      return (r.rotateAround(-2 * n, t), r.distanceSq(e) > i ? n : -n);
    }
    getRotateToAngle(t, n) {
      let r = e.subtract(this, n),
        i = e.subtract(t, n),
        a = e.angle(i, r);
      return Number.isNaN(a) ? 0 : a;
    }
    isInPolygon(e) {
      let t = !1;
      for (let n = 0, r = e.length - 1; n < e.length; r = n++) {
        let i = e[n],
          a = e[r];
        i.y > this.y != a.y > this.y &&
          this.x < ((a.x - i.x) * (this.y - i.y)) / (a.y - i.y) + i.x &&
          (t = !t);
      }
      return t;
    }
    length() {
      return Math.sqrt(this.x * this.x + this.y * this.y);
    }
    lengthSq() {
      return this.x * this.x + this.y * this.y;
    }
    normalize() {
      return (this.divide(this.length()), this);
    }
    normalized() {
      return e.divideScalar(this, this.length());
    }
    whichSide(e, t) {
      return (this.x - e.x) * (t.y - e.y) - (this.y - e.y) * (t.x - e.x);
    }
    sameSideAs(e, t, n) {
      let r = this.whichSide(e, t),
        i = n.whichSide(e, t);
      return (r < 0 && i < 0) || (r == 0 && i == 0) || (r > 0 && i > 0);
    }
    static add(t, n) {
      return new e(t.x + n.x, t.y + n.y);
    }
    static subtract(t, n) {
      return new e(t.x - n.x, t.y - n.y);
    }
    static multiply(t, n) {
      return new e(t.x * n.x, t.y * n.y);
    }
    static multiplyScalar(t, n) {
      return new e(t.x, t.y).multiplyScalar(n);
    }
    static midpoint(t, n) {
      return new e((t.x + n.x) / 2, (t.y + n.y) / 2);
    }
    static normals(t, n) {
      let r = e.subtract(n, t);
      return [new e(-r.y, r.x), new e(r.y, -r.x)];
    }
    static units(t, n) {
      let r = e.subtract(n, t);
      return [new e(-r.y, r.x).normalize(), new e(r.y, -r.x).normalize()];
    }
    static divide(t, n) {
      return new e(t.x / n.x, t.y / n.y);
    }
    static divideScalar(t, n) {
      return new e(t.x / n, t.y / n);
    }
    static dot(e, t) {
      return e.x * t.x + e.y * t.y;
    }
    static angle(t, n) {
      let r = e.dot(t, n);
      return Math.acos(r / (t.length() * n.length()));
    }
    static threePointangle(t, n, r) {
      let i = e.subtract(n, t),
        a = e.subtract(r, n),
        o = t.distance(n),
        s = n.distance(r);
      return Math.acos(e.dot(i, a) / (o * s));
    }
    static scalarProjection(t, n) {
      let r = n.normalized();
      return e.dot(t, r);
    }
    static averageDirection(t) {
      let n = new e(0, 0);
      for (let e = 0; e < t.length; e++) {
        let r = t[e];
        n.add(r);
      }
      return n.normalize();
    }
  },
  g = class e {
    constructor(
      e = new h(0, 0),
      t = new h(0, 0),
      n = null,
      r = null,
      i = !1,
      a = !1,
    ) {
      ((this.from = e),
        (this.to = t),
        (this.elementFrom = n),
        (this.elementTo = r),
        (this.chiralFrom = i),
        (this.chiralTo = a));
    }
    clone() {
      return new e(
        this.from.clone(),
        this.to.clone(),
        this.elementFrom,
        this.elementTo,
      );
    }
    getLength() {
      let e = this.to.x - this.from.x,
        t = this.to.y - this.from.y;
      return Math.sqrt(e * e + t * t);
    }
    getAngle() {
      return h.subtract(this.getRightVector(), this.getLeftVector()).angle();
    }
    getRightVector() {
      return this.from.x < this.to.x ? this.to : this.from;
    }
    getLeftVector() {
      return this.from.x < this.to.x ? this.from : this.to;
    }
    getRightElement() {
      return this.from.x < this.to.x ? this.elementTo : this.elementFrom;
    }
    getLeftElement() {
      return this.from.x < this.to.x ? this.elementFrom : this.elementTo;
    }
    getRightChiral() {
      return this.from.x < this.to.x ? this.chiralTo : this.chiralFrom;
    }
    getLeftChiral() {
      return this.from.x < this.to.x ? this.chiralFrom : this.chiralTo;
    }
    setRightVector(e, t) {
      return (
        this.from.x < this.to.x
          ? ((this.to.x = e), (this.to.y = t))
          : ((this.from.x = e), (this.from.y = t)),
        this
      );
    }
    setLeftVector(e, t) {
      return (
        this.from.x < this.to.x
          ? ((this.from.x = e), (this.from.y = t))
          : ((this.to.x = e), (this.to.y = t)),
        this
      );
    }
    rotateToXAxis() {
      let e = this.getLeftVector();
      return (this.setRightVector(e.x + this.getLength(), e.y), this);
    }
    rotate(e) {
      let t = this.getLeftVector(),
        n = this.getRightVector(),
        r = Math.sin(e),
        i = Math.cos(e),
        a = i * (n.x - t.x) - r * (n.y - t.y) + t.x,
        o = r * (n.x - t.x) - i * (n.y - t.y) + t.y;
      return (this.setRightVector(a, o), this);
    }
    shortenFrom(e) {
      let t = h.subtract(this.to, this.from);
      return (t.normalize(), t.multiplyScalar(e), this.from.add(t), this);
    }
    shortenTo(e) {
      let t = h.subtract(this.from, this.to);
      return (t.normalize(), t.multiplyScalar(e), this.to.add(t), this);
    }
    shortenRight(e) {
      return (
        this.from.x < this.to.x ? this.shortenTo(e) : this.shortenFrom(e),
        this
      );
    }
    shortenLeft(e) {
      return (
        this.from.x < this.to.x ? this.shortenFrom(e) : this.shortenTo(e),
        this
      );
    }
    shorten(e) {
      let t = h.subtract(this.from, this.to);
      return (
        t.normalize(),
        t.multiplyScalar(e / 2),
        this.to.add(t),
        this.from.subtract(t),
        this
      );
    }
  },
  _ = class e {
    static round(e, t) {
      if (t) {
        let n = 10 ** t;
        return Math.round(e * n) / n;
      }
      return Math.round(e);
    }
    static meanAngle(e) {
      let t = 0,
        n = 0;
      for (let r = 0; r < e.length; r++)
        ((t += Math.sin(e[r])), (n += Math.cos(e[r])));
      return Math.atan2(t / e.length, n / e.length);
    }
    static innerAngle(t) {
      return e.toRad(((t - 2) * 180) / t);
    }
    static polyCircumradius(e, t) {
      return e / (2 * Math.sin(Math.PI / t));
    }
    static apothem(e, t) {
      return e * Math.cos(Math.PI / t);
    }
    static apothemFromSideLength(t, n) {
      let r = e.polyCircumradius(t, n);
      return e.apothem(r, n);
    }
    static centralAngle(t) {
      return e.toRad(360 / t);
    }
    static toDeg(t) {
      return t * e.degFactor;
    }
    static toRad(t) {
      return t * e.radFactor;
    }
    static parityOfPermutation(e) {
      let t = new Uint8Array(e.length),
        n = 0,
        r = function (n, i = 0) {
          return t[n] === 1 ? i : (i++, (t[n] = 1), r(e[n], i));
        };
      for (let i = 0; i < e.length; i++) {
        if (t[i] === 1) continue;
        let e = r(i);
        n += 1 - (e % 2);
      }
      return n % 2 ? -1 : 1;
    }
    static get radFactor() {
      return Math.PI / 180;
    }
    static get degFactor() {
      return 180 / Math.PI;
    }
    static get twoPI() {
      return 2 * Math.PI;
    }
  },
  v = class e {
    constructor(e, t = 0, n = 0) {
      ((this.id = null),
        (this.value = e),
        (this.position = new h(t || 0, n || 0)),
        (this.previousPosition = new h(0, 0)),
        (this.parentVertexId = null),
        (this.children = []),
        (this.spanningTreeChildren = []),
        (this.edges = []),
        (this.positioned = !1),
        (this.angle = null),
        (this.dir = 1),
        (this.neighbourCount = 0),
        (this.neighbours = []),
        (this.neighbouringElements = []),
        (this.forcePositioned = !1));
    }
    setPosition(e, t) {
      ((this.position.x = e), (this.position.y = t));
    }
    setPositionFromVector(e) {
      ((this.position.x = e.x), (this.position.y = e.y));
    }
    addChild(e) {
      (this.children.push(e), this.neighbours.push(e), this.neighbourCount++);
    }
    addRingbondChild(e, t) {
      if ((this.children.push(e), this.value.bracket)) {
        let n = 1;
        (this.id === 0 && this.value.bracket.hcount === 0 && (n = 0),
          this.value.bracket.hcount === 1 && t === 0 && (n = 2),
          this.value.bracket.hcount === 1 &&
            t === 1 &&
            (n = this.neighbours.length < 3 ? 2 : 3),
          this.value.bracket.hcount === null && t === 0 && (n = 1),
          this.value.bracket.hcount === null &&
            t === 1 &&
            (n = this.neighbours.length < 3 ? 1 : 2),
          this.neighbours.splice(n, 0, e));
      } else this.neighbours.push(e);
      this.neighbourCount++;
    }
    setParentVertexId(e) {
      (this.neighbourCount++,
        (this.parentVertexId = e),
        this.neighbours.push(e));
    }
    isTerminal() {
      return this.value.hasAttachedPseudoElements
        ? !0
        : (this.parentVertexId === null && this.children.length < 2) ||
            this.children.length === 0;
    }
    clone() {
      let t = new e(this.value, this.position.x, this.position.y);
      return (
        (t.id = this.id),
        (t.previousPosition = new h(
          this.previousPosition.x,
          this.previousPosition.y,
        )),
        (t.parentVertexId = this.parentVertexId),
        (t.children = f.clone(this.children)),
        (t.spanningTreeChildren = f.clone(this.spanningTreeChildren)),
        (t.edges = f.clone(this.edges)),
        (t.positioned = this.positioned),
        (t.angle = this.angle),
        (t.forcePositioned = this.forcePositioned),
        t
      );
    }
    equals(e) {
      return this.id === e.id;
    }
    getAngle(e = null, t = !1) {
      let n = null;
      return (
        (n = e
          ? h.subtract(this.position, e)
          : h.subtract(this.position, this.previousPosition)),
        t ? _.toDeg(n.angle()) : n.angle()
      );
    }
    getTextDirection(e, t = !1) {
      let n = this.getDrawnNeighbours(e),
        r = [];
      if (e.length === 1) return `right`;
      for (let t = 0; t < n.length; t++)
        r.push(this.getAngle(e[n[t]].position));
      let i = _.meanAngle(r);
      if (this.isTerminal() || t)
        (Math.round(i * 100) / 100 == 1.57 && (i -= 0.2),
          (i = Math.round(Math.round(i / Math.PI) * Math.PI)));
      else {
        let e = Math.PI / 2;
        i = Math.round(Math.round(i / e) * e);
      }
      return i === 2
        ? `down`
        : i === -2
          ? `up`
          : i === 0
            ? `right`
            : i === 3 || i === -3
              ? `left`
              : `down`;
    }
    getNeighbours(e = null) {
      if (e === null) return this.neighbours.slice();
      let t = [];
      for (let n = 0; n < this.neighbours.length; n++)
        this.neighbours[n] !== e && t.push(this.neighbours[n]);
      return t;
    }
    getDrawnNeighbours(e) {
      let t = [];
      for (let n = 0; n < this.neighbours.length; n++)
        e[this.neighbours[n]].value.isDrawn && t.push(this.neighbours[n]);
      return t;
    }
    getNeighbourCount() {
      return this.neighbourCount;
    }
    getSpanningTreeNeighbours(e = null) {
      let t = [];
      for (let n = 0; n < this.spanningTreeChildren.length; n++)
        (e === void 0 || e != this.spanningTreeChildren[n]) &&
          t.push(this.spanningTreeChildren[n]);
      return (
        this.parentVertexId != null &&
          (e === void 0 || e != this.parentVertexId) &&
          t.push(this.parentVertexId),
        t
      );
    }
    getNextInRing(e, t, n) {
      let r = this.getNeighbours();
      for (let i = 0; i < r.length; i++)
        if (f.contains(e[r[i]].value.rings, { value: t }) && r[i] != n)
          return r[i];
      return null;
    }
  },
  y = class {
    constructor(e, t) {
      ((this.id = null),
        (this.firstRingId = e.id),
        (this.secondRingId = t.id),
        (this.vertices = new Set()),
        (this.isForcedBridge = !1));
      for (let n = 0; n < e.members.length; n++) {
        let r = e.members[n];
        for (let e = 0; e < t.members.length; e++)
          r === t.members[e] && this.addVertex(r);
      }
    }
    addVertex(e) {
      this.vertices.add(e);
    }
    updateOther(e, t) {
      this.firstRingId === t ? (this.secondRingId = e) : (this.firstRingId = e);
    }
    containsRing(e) {
      return this.firstRingId === e || this.secondRingId === e;
    }
    isBridge(e) {
      if (this.isForcedBridge || this.vertices.size > 2) return !0;
      if (this.vertices.size === 2) {
        let [t, n] = [...this.vertices],
          r = new Set(e[n].neighbours);
        for (let i of e[t].neighbours)
          if (i !== n && r.has(i)) {
            let t = e[i].value.rings;
            if (t.includes(this.firstRingId) || t.includes(this.secondRingId))
              return !0;
          }
      }
      return !1;
    }
    static isBridge(e, t, n, r) {
      let i = null;
      for (let a = 0; a < e.length; a++)
        if (
          ((i = e[a]),
          (i.firstRingId === n && i.secondRingId === r) ||
            (i.firstRingId === r && i.secondRingId === n))
        )
          return i.isBridge(t);
      return !1;
    }
    static getNeighbours(e, t) {
      let n = [];
      for (let r = 0; r < e.length; r++) {
        let i = e[r];
        i.firstRingId === t
          ? n.push(i.secondRingId)
          : i.secondRingId === t && n.push(i.firstRingId);
      }
      return n;
    }
    static getVertices(e, t, n) {
      for (let r = 0; r < e.length; r++) {
        let i = e[r];
        if (
          (i.firstRingId === t && i.secondRingId === n) ||
          (i.firstRingId === n && i.secondRingId === t)
        )
          return [...i.vertices];
      }
    }
  },
  b = class e {
    constructor(e) {
      ((this.id = null),
        (this.members = e),
        (this.edges = []),
        (this.insiders = []),
        (this.neighbours = []),
        (this.positioned = !1),
        (this.center = new h(0, 0)),
        (this.rings = []),
        (this.isBridged = !1),
        (this.isPartOfBridged = !1),
        (this.isSpiro = !1),
        (this.isFused = !1),
        (this.centralAngle = 0),
        (this.canFlip = !0));
    }
    clone() {
      let t = new e(this.members);
      return (
        (t.id = this.id),
        (t.insiders = f.clone(this.insiders)),
        (t.neighbours = f.clone(this.neighbours)),
        (t.positioned = this.positioned),
        (t.center = this.center.clone()),
        (t.rings = f.clone(this.rings)),
        (t.isBridged = this.isBridged),
        (t.isPartOfBridged = this.isPartOfBridged),
        (t.isSpiro = this.isSpiro),
        (t.isFused = this.isFused),
        (t.centralAngle = this.centralAngle),
        (t.canFlip = this.canFlip),
        t
      );
    }
    getSize() {
      return this.members.length;
    }
    getPolygon(e) {
      let t = [];
      for (let n = 0; n < this.members.length; n++)
        t.push(e[this.members[n]].position);
      return t;
    }
    getAngle() {
      return Math.PI - this.centralAngle;
    }
    eachMember(e, t, n, r) {
      n = n || n === 0 ? n : this.members[0];
      let i = n,
        a = 0;
      for (; i != null && a < 100;) {
        let o = i;
        (t(o),
          (i = e[i].getNextInRing(e, this.id, r)),
          (r = o),
          i == n && (i = null),
          a++);
      }
    }
    getOrderedNeighbours(e) {
      let t = Array(this.neighbours.length);
      for (let n = 0; n < this.neighbours.length; n++) {
        let r = y.getVertices(e, this.id, this.neighbours[n]);
        t[n] = { n: r.length, neighbour: this.neighbours[n] };
      }
      return (t.sort((e, t) => t.n - e.n), t);
    }
    isBenzeneLike(e) {
      let t = this.getDoubleBondCount(e),
        n = this.members.length;
      return (t === 3 && n === 6) || (t === 2 && n === 5);
    }
    getDoubleBondCount(e) {
      let t = 0;
      for (let n = 0; n < this.members.length; n++) {
        let r = e[this.members[n]].value;
        (r.bondType === `=` || r.branchBond === `=`) && t++;
      }
      return t;
    }
    contains(e) {
      for (let t = 0; t < this.members.length; t++)
        if (this.members[t] == e) return !0;
      return !1;
    }
  },
  x = class {
    constructor(e, t) {
      ((this.colors = e), (this.theme = this.colors[t]));
    }
    getColor(e) {
      return e && ((e = e.toUpperCase()), e in this.theme)
        ? this.theme[e]
        : this.theme.C;
    }
    setTheme(e) {
      e in this.colors && (this.theme = this.colors[e]);
    }
  };
function S(e) {
  return e ? (e === 1 ? `+` : e === -1 ? `-` : e > 0 ? e + `+` : e + `-`) : ``;
}
var C = class {
    constructor(e, t, n) {
      let r = null;
      if (
        ((r =
          e instanceof String
            ? document.getElementById(e.valueOf())
            : typeof e == `string`
              ? document.getElementById(e)
              : e),
        r instanceof HTMLCanvasElement)
      )
        this.canvas = r;
      else
        throw Error(`First argument was not a canvas or the ID of a canvas.`);
      ((this.ctx = this.canvas.getContext(`2d`)),
        (this.themeManager = t),
        (this.opts = n),
        (this.drawingWidth = 0),
        (this.drawingHeight = 0),
        (this.offsetX = 0),
        (this.offsetY = 0),
        (this.fontLarge =
          this.opts.fontSizeLarge + `pt Helvetica, Arial, sans-serif`),
        (this.fontSmall =
          this.opts.fontSizeSmall + `pt Helvetica, Arial, sans-serif`),
        this.updateSize(this.opts.width, this.opts.height),
        (this.ctx.font = this.fontLarge),
        (this.hydrogenWidth = this.ctx.measureText(`H`).width),
        (this.halfHydrogenWidth = this.hydrogenWidth / 2),
        (this.halfBondThickness = this.opts.bondThickness / 2));
    }
    updateSize(e, t) {
      ((this.ratio = window.devicePixelRatio || 1),
        this.ratio === 1
          ? ((this.canvas.width = e * this.ratio),
            (this.canvas.height = t * this.ratio))
          : ((this.canvas.width = e * this.ratio),
            (this.canvas.height = t * this.ratio),
            (this.canvas.style.width = e + `px`),
            (this.canvas.style.height = t + `px`),
            this.ctx.setTransform(this.ratio, 0, 0, this.ratio, 0, 0)));
    }
    setTheme(e) {
      this.colors = e;
    }
    scale(e) {
      let t = -Number.MAX_VALUE,
        n = -Number.MAX_VALUE,
        r = Number.MAX_VALUE,
        i = Number.MAX_VALUE;
      for (let a = 0; a < e.length; a++) {
        if (!e[a].value.isDrawn) continue;
        let o = e[a].position;
        (t < o.x && (t = o.x),
          n < o.y && (n = o.y),
          r > o.x && (r = o.x),
          i > o.y && (i = o.y));
      }
      let a = this.opts.padding;
      ((t += a),
        (n += a),
        (r -= a),
        (i -= a),
        (this.drawingWidth = t - r),
        (this.drawingHeight = n - i));
      let o = this.canvas.offsetWidth / this.drawingWidth,
        s = this.canvas.offsetHeight / this.drawingHeight,
        c = o < s ? o : s;
      (this.ctx.scale(c, c),
        (this.offsetX = -r),
        (this.offsetY = -i),
        o < s
          ? (this.offsetY +=
              this.canvas.offsetHeight / (2 * c) - this.drawingHeight / 2)
          : (this.offsetX +=
              this.canvas.offsetWidth / (2 * c) - this.drawingWidth / 2));
    }
    reset() {
      this.ctx.setTransform(1, 0, 0, 1, 0, 0);
    }
    getColor(e) {
      return (
        (e = e.toUpperCase()),
        e in this.colors ? this.colors[e] : this.colors.C
      );
    }
    drawCircle(e, t, n, r, i = !0, a = !1, o = ``) {
      let s = this.ctx,
        c = this.offsetX,
        l = this.offsetY;
      (s.save(),
        (s.lineWidth = 1.5),
        s.beginPath(),
        s.arc(e + c, t + l, n, 0, _.twoPI, !0),
        s.closePath(),
        a
          ? (i
              ? ((s.fillStyle = `#f00`), s.fill())
              : ((s.strokeStyle = `#f00`), s.stroke()),
            this.drawDebugText(e, t, o))
          : i
            ? ((s.fillStyle = r), s.fill())
            : ((s.strokeStyle = r), s.stroke()),
        s.restore());
    }
    drawLine(e, t = !1, n = 1) {
      let r = this.ctx,
        i = this.offsetX,
        a = this.offsetY,
        o = e.clone().shorten(4),
        s = o.getLeftVector().clone(),
        c = o.getRightVector().clone();
      ((s.x += i),
        (s.y += a),
        (c.x += i),
        (c.y += a),
        t ||
          (r.save(),
          (r.globalCompositeOperation = `destination-out`),
          r.beginPath(),
          r.moveTo(s.x, s.y),
          r.lineTo(c.x, c.y),
          (r.lineCap = `round`),
          (r.lineWidth = this.opts.bondThickness + 1.2),
          (r.strokeStyle = this.themeManager.getColor(`BACKGROUND`)),
          r.stroke(),
          (r.globalCompositeOperation = `source-over`),
          r.restore()),
        (s = e.getLeftVector().clone()),
        (c = e.getRightVector().clone()),
        (s.x += i),
        (s.y += a),
        (c.x += i),
        (c.y += a),
        r.save(),
        r.beginPath(),
        r.moveTo(s.x, s.y),
        r.lineTo(c.x, c.y),
        (r.lineCap = `round`),
        (r.lineWidth = this.opts.bondThickness));
      let l = this.ctx.createLinearGradient(s.x, s.y, c.x, c.y);
      (l.addColorStop(
        0.4,
        this.themeManager.getColor(e.getLeftElement()) ||
          this.themeManager.getColor(`C`),
      ),
        l.addColorStop(
          0.6,
          this.themeManager.getColor(e.getRightElement()) ||
            this.themeManager.getColor(`C`),
        ),
        t &&
          (r.setLineDash([1, 1.5]),
          (r.lineWidth = this.opts.bondThickness / 1.5)),
        n < 1 && (r.globalAlpha = n),
        (r.strokeStyle = l),
        r.stroke(),
        r.restore());
    }
    drawWedge(e, t = 1) {
      if (isNaN(e.from.x) || isNaN(e.from.y) || isNaN(e.to.x) || isNaN(e.to.y))
        return;
      let n = this.ctx,
        r = this.offsetX,
        i = this.offsetY,
        a = e.clone().shorten(5),
        o = a.getLeftVector().clone(),
        s = a.getRightVector().clone();
      ((o.x += r),
        (o.y += i),
        (s.x += r),
        (s.y += i),
        (o = e.getLeftVector().clone()),
        (s = e.getRightVector().clone()),
        (o.x += r),
        (o.y += i),
        (s.x += r),
        (s.y += i),
        n.save());
      let c = h.normals(o, s);
      (c[0].normalize(), c[1].normalize());
      let l = e.getRightChiral(),
        u = o,
        d = s;
      l && ((u = s), (d = o));
      let f = h.add(u, h.multiplyScalar(c[0], this.halfBondThickness)),
        p = h.add(d, h.multiplyScalar(c[0], 1.5 + this.halfBondThickness)),
        m = h.add(d, h.multiplyScalar(c[1], 1.5 + this.halfBondThickness)),
        g = h.add(u, h.multiplyScalar(c[1], this.halfBondThickness));
      (n.beginPath(),
        n.moveTo(f.x, f.y),
        n.lineTo(p.x, p.y),
        n.lineTo(m.x, m.y),
        n.lineTo(g.x, g.y));
      let _ = this.ctx.createRadialGradient(
        s.x,
        s.y,
        this.opts.bondLength,
        s.x,
        s.y,
        0,
      );
      (_.addColorStop(
        0.4,
        this.themeManager.getColor(e.getLeftElement()) ||
          this.themeManager.getColor(`C`),
      ),
        _.addColorStop(
          0.6,
          this.themeManager.getColor(e.getRightElement()) ||
            this.themeManager.getColor(`C`),
        ),
        (n.fillStyle = _),
        n.fill(),
        n.restore());
    }
    drawDashedWedge(e) {
      if (isNaN(e.from.x) || isNaN(e.from.y) || isNaN(e.to.x) || isNaN(e.to.y))
        return;
      let t = this.ctx,
        n = this.offsetX,
        r = this.offsetY,
        i = e.getLeftVector().clone(),
        a = e.getRightVector().clone();
      ((i.x += n), (i.y += r), (a.x += n), (a.y += r), t.save());
      let o = h.normals(i, a);
      (o[0].normalize(), o[1].normalize());
      let s = e.getRightChiral(),
        c,
        l,
        u,
        d,
        f = e.clone();
      (s
        ? ((c = a),
          (l = i),
          f.shortenRight(1),
          (u = f.getRightVector().clone()),
          (d = f.getLeftVector().clone()))
        : ((c = i),
          (l = a),
          f.shortenLeft(1),
          (u = f.getLeftVector().clone()),
          (d = f.getRightVector().clone())),
        (u.x += n),
        (u.y += r),
        (d.x += n),
        (d.y += r));
      let p = h.subtract(l, c).normalize();
      ((t.strokeStyle = this.themeManager.getColor(`C`)),
        (t.lineCap = `round`),
        (t.lineWidth = this.opts.bondThickness),
        t.beginPath());
      let m = e.getLength(),
        g = 1.25 / (m / (this.opts.bondThickness * 3)),
        _ = !1;
      for (let n = 0; n < 1; n += g) {
        let r = h.multiplyScalar(p, n * m),
          i = h.add(c, r),
          a = 1.5 * n,
          s = h.multiplyScalar(o[0], a);
        (!_ &&
          n > 0.5 &&
          (t.stroke(),
          t.beginPath(),
          (t.strokeStyle =
            this.themeManager.getColor(e.getRightElement()) ||
            this.themeManager.getColor(`C`)),
          (_ = !0)),
          i.subtract(s),
          t.moveTo(i.x, i.y),
          i.add(h.multiplyScalar(s, 2)),
          t.lineTo(i.x, i.y));
      }
      (t.stroke(), t.restore());
    }
    drawDebugText(e, t, n) {
      let r = this.ctx;
      (r.save(),
        (r.font = `5px Droid Sans, sans-serif`),
        (r.textAlign = `start`),
        (r.textBaseline = `top`),
        (r.fillStyle = `#ff0000`),
        r.fillText(n, e + this.offsetX, t + this.offsetY),
        r.restore());
    }
    drawBall(e, t, n) {
      let r = this.ctx;
      (r.save(),
        r.beginPath(),
        r.arc(
          e + this.offsetX,
          t + this.offsetY,
          this.opts.bondLength / 4.5,
          0,
          _.twoPI,
          !1,
        ),
        (r.fillStyle = this.themeManager.getColor(n)),
        r.fill(),
        r.restore());
    }
    drawPoint(e, t, n) {
      let r = this.ctx,
        i = this.offsetX,
        a = this.offsetY;
      (r.save(),
        (r.globalCompositeOperation = `destination-out`),
        r.beginPath(),
        r.arc(e + i, t + a, 1.5, 0, _.twoPI, !0),
        r.closePath(),
        r.fill(),
        (r.globalCompositeOperation = `source-over`),
        r.beginPath(),
        r.arc(e + this.offsetX, t + this.offsetY, 0.75, 0, _.twoPI, !1),
        (r.fillStyle = this.themeManager.getColor(n)),
        r.fill(),
        r.restore());
    }
    drawText(e, t, n, r, i, a, o, s, c, l = {}) {
      let u = this.ctx,
        d = this.offsetX,
        f = this.offsetY;
      (u.save(), (u.textAlign = `start`), (u.textBaseline = `alphabetic`));
      let p = ``,
        m = 0;
      o &&
        ((p = S(o)), (u.font = this.fontSmall), (m = u.measureText(p).width));
      let h = `0`,
        g = 0;
      (s > 0 &&
        ((h = s.toString()),
        (u.font = this.fontSmall),
        (g = u.measureText(h).width)),
        o === 1 &&
          n === `N` &&
          `0O` in l &&
          `0O-1` in l &&
          ((l = {
            "0O": {
              element: `O`,
              count: 2,
              hydrogenCount: 0,
              previousElement: `C`,
              charge: ``,
            },
          }),
          (o = 0)),
        (u.font = this.fontLarge),
        (u.fillStyle = this.themeManager.getColor(`BACKGROUND`)));
      let v = u.measureText(n),
        y =
          v.width > this.opts.fontSizeLarge ? v.width : this.opts.fontSizeLarge;
      ((y /= 1.5),
        (u.globalCompositeOperation = `destination-out`),
        u.beginPath(),
        u.arc(e + d, t + f, y, 0, _.twoPI, !0),
        u.closePath(),
        u.fill(),
        (u.globalCompositeOperation = `source-over`));
      let b = -v.width / 2,
        x = -v.width / 2;
      ((u.fillStyle = this.themeManager.getColor(n)),
        u.fillText(n, e + d + b, t + this.opts.halfFontSizeLarge + f),
        (b += v.width),
        o &&
          ((u.font = this.fontSmall),
          u.fillText(p, e + d + b, t - this.opts.fifthFontSizeSmall + f),
          (b += m)),
        s > 0 &&
          ((u.font = this.fontSmall),
          u.fillText(h, e + d + x - g, t - this.opts.fifthFontSizeSmall + f),
          (x -= g)),
        (u.font = this.fontLarge));
      let C = 0,
        w = 0;
      if (r === 1) {
        let n = e + d,
          r = t + f + this.opts.halfFontSizeLarge;
        ((C = this.hydrogenWidth),
          (x -= C),
          i === `left`
            ? (n += x)
            : i === `right` || (i === `up` && a) || (i === `down` && a)
              ? (n += b)
              : i === `up` && !a
                ? ((r -=
                    this.opts.fontSizeLarge + this.opts.quarterFontSizeLarge),
                  (n -= this.halfHydrogenWidth))
                : i === `down` &&
                  !a &&
                  ((r +=
                    this.opts.fontSizeLarge + this.opts.quarterFontSizeLarge),
                  (n -= this.halfHydrogenWidth)),
          u.fillText(`H`, n, r),
          (b += C));
      } else if (r > 1) {
        let n = e + d,
          o = t + f + this.opts.halfFontSizeLarge;
        ((C = this.hydrogenWidth),
          (u.font = this.fontSmall),
          (w = u.measureText(r.toString()).width),
          (x -= C + w),
          i === `left`
            ? (n += x)
            : i === `right` || (i === `up` && a) || (i === `down` && a)
              ? (n += b)
              : i === `up` && !a
                ? ((o -=
                    this.opts.fontSizeLarge + this.opts.quarterFontSizeLarge),
                  (n -= this.halfHydrogenWidth))
                : i === `down` &&
                  !a &&
                  ((o +=
                    this.opts.fontSizeLarge + this.opts.quarterFontSizeLarge),
                  (n -= this.halfHydrogenWidth)),
          (u.font = this.fontLarge),
          u.fillText(`H`, n, o),
          (u.font = this.fontSmall),
          u.fillText(
            r.toString(),
            n + this.halfHydrogenWidth + w,
            o + this.opts.fifthFontSizeSmall,
          ),
          (b += C + this.halfHydrogenWidth + w));
      }
      for (let n of Object.keys(l)) {
        let r = 0,
          a = 0,
          o = l[n].element,
          s = l[n].count,
          c = l[n].hydrogenCount,
          p = l[n].charge;
        ((u.font = this.fontLarge),
          s > 1 &&
            c > 0 &&
            ((r = u.measureText(`(`).width), (a = u.measureText(`)`).width)));
        let m = u.measureText(o).width,
          h = 0,
          g = ``,
          _ = 0;
        ((C = 0),
          c > 0 && (C = this.hydrogenWidth),
          (u.font = this.fontSmall),
          s > 1 && (h = u.measureText(s).width),
          p !== 0 && ((g = S(p)), (_ = u.measureText(g).width)),
          (w = 0),
          c > 1 && (w = u.measureText(c).width),
          (u.font = this.fontLarge));
        let v = e + d,
          y = t + f + this.opts.halfFontSizeLarge;
        ((u.fillStyle = this.themeManager.getColor(o)),
          s > 0 && (x -= h),
          s > 1 &&
            c > 0 &&
            (i === `left`
              ? ((x -= a), u.fillText(`)`, v + x, y))
              : (u.fillText(`(`, v + b, y), (b += r))),
          i === `left`
            ? ((x -= m), u.fillText(o, v + x, y))
            : (u.fillText(o, v + b, y), (b += m)),
          c > 0 &&
            (i === `left`
              ? ((x -= C + w),
                u.fillText(`H`, v + x, y),
                c > 1 &&
                  ((u.font = this.fontSmall),
                  u.fillText(c, v + x + C, y + this.opts.fifthFontSizeSmall)))
              : (u.fillText(`H`, v + b, y),
                (b += C),
                c > 1 &&
                  ((u.font = this.fontSmall),
                  u.fillText(c, v + b, y + this.opts.fifthFontSizeSmall),
                  (b += w)))),
          (u.font = this.fontLarge),
          s > 1 &&
            c > 0 &&
            (i === `left`
              ? ((x -= r), u.fillText(`(`, v + x, y))
              : (u.fillText(`)`, v + b, y), (b += a))),
          (u.font = this.fontSmall),
          s > 1 &&
            (i === `left`
              ? u.fillText(
                  s,
                  v + x + r + a + C + w + m,
                  y + this.opts.fifthFontSizeSmall,
                )
              : (u.fillText(s, v + b, y + this.opts.fifthFontSizeSmall),
                (b += h))),
          p !== 0 &&
            (i === `left`
              ? u.fillText(
                  g,
                  v + x + r + a + C + w + m,
                  t - this.opts.fifthFontSizeSmall + f,
                )
              : (u.fillText(g, v + b, t - this.opts.fifthFontSizeSmall + f),
                (b += _))));
      }
      u.restore();
    }
    drawDebugPoint(e, t, n = ``, r = `#f00`) {
      this.drawCircle(e, t, 2, r, !0, !0, n);
    }
    drawAromaticityRing(e) {
      let t = this.ctx,
        n = _.apothemFromSideLength(this.opts.bondLength, e.getSize());
      (t.save(),
        (t.strokeStyle = this.themeManager.getColor(`C`)),
        (t.lineWidth = this.opts.bondThickness),
        t.beginPath(),
        t.arc(
          e.center.x + this.offsetX,
          e.center.y + this.offsetY,
          n - this.opts.bondSpacing,
          0,
          Math.PI * 2,
          !0,
        ),
        t.closePath(),
        t.stroke(),
        t.restore());
    }
    clear() {
      this.ctx.clearRect(
        0,
        0,
        this.canvas.offsetWidth,
        this.canvas.offsetHeight,
      );
    }
  },
  w = class e {
    static build(t, n) {
      let r = new e(t, n, null, new Map());
      return (r.findChildren(), r.sortChildren(), r);
    }
    static compareAtoms(e, t) {
      if (e.atomicNumber !== t.atomicNumber)
        return [t.atomicNumber - e.atomicNumber, !1];
      if (e.cloneDepth && t.cloneDepth) {
        let n = e.cloneDepth - t.cloneDepth;
        if (n !== 0) return [n, !1];
      }
      return e.atomicWeight === t.atomicWeight
        ? [0, e.stereocenter || t.stereocenter]
        : [t.atomicWeight - e.atomicWeight, !1];
    }
    static compareTrees(e, t) {
      let n = this.compareAtoms(e, t);
      if (n[0] !== 0) return n;
      let r = [e],
        i = [t],
        a = !1;
      for (; r.length !== 0;) {
        for (let e = 0; e < r.length; ++e) {
          let t = r[e].findChildren(),
            o = i[e].findChildren(),
            s = Math.max(t.length, o.length);
          for (let e = 0; e < s; ++e) {
            if (e >= t.length) return [1, !1];
            if (e >= o.length) return [-1, !1];
            if (((n = this.compareAtoms(t[e], o[e])), n[0] !== 0)) return n;
            a ||= n[1];
          }
        }
        ((r = r.flatMap((e) => e.sortChildren())),
          (i = i.flatMap((e) => e.sortChildren())));
      }
      return [0, a];
    }
    constructor(e, t, n, r) {
      ((this.graph = e),
        (this.parent = n),
        t instanceof v
          ? ((this.vertex = t),
            (this.atomicNumber = t.value.getAtomicNumber()),
            (this.atomicWeight = 0),
            (this.stereocenter = !!(
              t.value.bracket && t.value.bracket.chirality
            )))
          : ((this.vertex = null),
            (this.atomicNumber = t),
            (this.atomicWeight = 0),
            (this.stereocenter = !1)),
        r instanceof Map
          ? ((this.cloneDepth = void 0),
            (this.visited = r),
            (this.children = void 0),
            (this.sorted = !1))
          : ((this.cloneDepth = r),
            (this.visited = null),
            (this.children = []),
            (this.sorted = !0),
            (this.stereocenter = !1)));
    }
    findChildren() {
      if (this.children === void 0) {
        this.children = [];
        let t = this.visited.size + 1,
          n = new Map(this.visited);
        (n.set(this.vertex, t), (this.visited = null));
        for (let r of this.vertex.neighbours) {
          let i = this.graph.getEdge(this.vertex.id, r).weight;
          if (!i) continue;
          let a = this.graph.vertices[r],
            o = n.get(a);
          if (Object.is(a, this.parent)) --i;
          else if (o === void 0) {
            let r = new e(this.graph, a, this.vertex, n);
            (this.children.push(r), (o = t + 1), --i);
          }
          if (i > 0) {
            let t = new e(this.graph, a, this.vertex, o);
            for (; i-- > 0;) this.children.push(t);
          }
        }
        if (this.vertex.value.isPartOfAromaticRing) {
          let t = new e(this.graph, 6, this.vertex);
          this.children.push(t);
        }
        let r = this.vertex.value.countImplicitHydrogens();
        if (r > 0) {
          let t = new e(this.graph, 1, this.vertex);
          for (; r-- > 0;) this.children.push(t);
        }
        this.children.sort((t, n) => e.compareAtoms(t, n)[0]);
      }
      return this.children;
    }
    size() {
      return this.children === void 0
        ? 1
        : this.children.reduce((e, t) => e + t.size(), 1);
    }
    sortChildren() {
      return (
        (this.sorted ||=
          (this.children.sort((t, n) => {
            let [r, i] = e.compareTrees(t, n);
            return i && this.stereocenter
              ? this.vertex.value.bracket.chirality === `@`
                ? n.vertex.id - t.vertex.id
                : t.vertex.id - n.vertex.id
              : r;
          }),
          !0)),
        this.children
      );
    }
  },
  T = class {
    static getOrderArray(e, t) {
      let n = w.build(e, t);
      for (let e = 1; e < n.children.length; ++e) {
        let t = w.compareTrees(n.children[e - 1], n.children[e]);
        if (t[0] === 0 && !t[1]) return;
      }
      let r = t.neighbours.length,
        i = Array(r);
      for (let e = 0; e < r; ++e) {
        let r = n.children[e].vertex.id;
        i[e] = t.neighbours.findIndex((e) => e === r);
      }
      return i;
    }
  },
  E = class e {
    constructor(e, t, n = 1) {
      ((this.id = null),
        (this.sourceId = e),
        (this.targetId = t),
        (this.weight = n),
        (this.bondType = `-`),
        (this.isPartOfAromaticRing = !1),
        (this.center = !1),
        (this.wedge = ``));
    }
    setBondType(t) {
      ((this.bondType = t), (this.weight = e.bonds[t]));
    }
    static get bonds() {
      return { ".": 0, "-": 1, "/": 1, "\\": 1, "=": 2, "#": 3, $: 4 };
    }
  },
  D = class e {
    constructor(e, t = !1) {
      ((this.vertices = []),
        (this.edges = []),
        (this.atomIdxToVertexId = []),
        (this.vertexIdsToEdgeId = {}),
        (this.isomeric = t),
        (this._atomIdx = 0),
        (this._time = 0),
        this._init(e));
    }
    _init(e, t = 0, n = null, r = !1) {
      let i = e.atom.element ? e.atom.element : e.atom,
        a = new m(i, e.bond);
      ((i !== `H` || (!e.hasNext && n === null)) &&
        ((a.idx = this._atomIdx), this._atomIdx++),
        (a.branchBond = e.branchBond),
        (a.ringbonds = e.ringbonds),
        (a.bracket = e.atom.element ? e.atom : null),
        (a.class = e.atom.class));
      let o = new v(a),
        s = this.vertices[n];
      if (
        (this.addVertex(o),
        a.idx !== null && this.atomIdxToVertexId.push(o.id),
        n !== null)
      ) {
        (o.setParentVertexId(n),
          o.value.addNeighbouringElement(s.value.element),
          s.addChild(o.id),
          s.value.addNeighbouringElement(a.element),
          s.spanningTreeChildren.push(o.id));
        let e = new E(n, o.id, 1);
        (r
          ? e.setBondType(o.value.branchBond || `-`)
          : e.setBondType(s.value.bondType || `-`),
          this.addEdge(e));
      }
      let c = e.ringbondCount + 1;
      a.bracket && (c += a.bracket.hcount);
      let l = 0;
      if (a.bracket && a.bracket.chirality) {
        ((a.isStereoCenter = !0), (l = a.bracket.hcount));
        for (let e = 0; e < l; e++)
          this._init(
            {
              atom: `H`,
              isBracket: `false`,
              branches: [],
              branchCount: 0,
              ringbonds: [],
              ringbondCount: !1,
              next: null,
              hasNext: !1,
              bond: `-`,
            },
            e,
            o.id,
            !0,
          );
      }
      for (let t = 0; t < e.branchCount; t++)
        this._init(e.branches[t], t + c, o.id, !0);
      e.hasNext && this._init(e.next, e.branchCount + c, o.id);
    }
    clear() {
      ((this.vertices = []), (this.edges = []), (this.vertexIdsToEdgeId = {}));
    }
    addVertex(e) {
      return ((e.id = this.vertices.length), this.vertices.push(e), e.id);
    }
    addEdge(e) {
      let t = this.vertices[e.sourceId],
        n = this.vertices[e.targetId];
      return (
        (e.id = this.edges.length),
        this.edges.push(e),
        (this.vertexIdsToEdgeId[e.sourceId + `_` + e.targetId] = e.id),
        (this.vertexIdsToEdgeId[e.targetId + `_` + e.sourceId] = e.id),
        (e.isPartOfAromaticRing =
          t.value.isPartOfAromaticRing && n.value.isPartOfAromaticRing),
        (t.value.bondCount += e.weight),
        (n.value.bondCount += e.weight),
        t.edges.push(e.id),
        n.edges.push(e.id),
        e.id
      );
    }
    getEdge(e, t) {
      let n = this.vertexIdsToEdgeId[e + `_` + t];
      return n === void 0 ? null : this.edges[n];
    }
    getEdges(e) {
      let t = [],
        n = this.vertices[e];
      for (let r = 0; r < n.neighbours.length; r++)
        t.push(this.vertexIdsToEdgeId[e + `_` + n.neighbours[r]]);
      return t;
    }
    hasEdge(e, t) {
      return this.vertexIdsToEdgeId[e + `_` + t] !== void 0;
    }
    getVertexList() {
      let e = [this.vertices.length];
      for (let t = 0; t < this.vertices.length; t++) e[t] = this.vertices[t].id;
      return e;
    }
    getEdgeList() {
      let e = Array(this.edges.length);
      for (let t = 0; t < this.edges.length; t++)
        e[t] = [this.edges[t].sourceId, this.edges[t].targetId];
      return e;
    }
    getAdjacencyMatrix() {
      let e = this.vertices.length,
        t = Array(e);
      for (let n = 0; n < e; n++) ((t[n] = Array(e)), t[n].fill(0));
      for (let e = 0; e < this.edges.length; e++) {
        let n = this.edges[e];
        ((t[n.sourceId][n.targetId] = 1), (t[n.targetId][n.sourceId] = 1));
      }
      return t;
    }
    getComponentsAdjacencyMatrix() {
      let e = this.vertices.length,
        t = Array(e),
        n = this.getBridges();
      for (let n = 0; n < e; n++) ((t[n] = Array(e)), t[n].fill(0));
      for (let e = 0; e < this.edges.length; e++) {
        let n = this.edges[e];
        ((t[n.sourceId][n.targetId] = 1), (t[n.targetId][n.sourceId] = 1));
      }
      for (let e = 0; e < n.length; e++)
        ((t[n[e][0]][n[e][1]] = 0), (t[n[e][1]][n[e][0]] = 0));
      return t;
    }
    getSubgraphAdjacencyMatrix(e) {
      let t = e.length,
        n = Array(t);
      for (let r = 0; r < t; r++) {
        ((n[r] = Array(t)), n[r].fill(0));
        for (let i = 0; i < t; i++)
          r !== i && this.hasEdge(e[r], e[i]) && (n[r][i] = 1);
      }
      return n;
    }
    getDistanceMatrix() {
      let e = this.vertices.length,
        t = this.getAdjacencyMatrix(),
        n = Array(e);
      for (let t = 0; t < e; t++) ((n[t] = Array(e)), n[t].fill(1 / 0));
      for (let r = 0; r < e; r++)
        for (let i = 0; i < e; i++) t[r][i] === 1 && (n[r][i] = 1);
      for (let t = 0; t < e; t++)
        for (let r = 0; r < e; r++)
          for (let i = 0; i < e; i++)
            n[r][i] > n[r][t] + n[t][i] && (n[r][i] = n[r][t] + n[t][i]);
      return n;
    }
    getSubgraphDistanceMatrix(e) {
      let t = e.length,
        n = this.getSubgraphAdjacencyMatrix(e),
        r = Array(t);
      for (let e = 0; e < t; e++)
        ((r[e] = Array(t)), r[e].fill(1 / 0), (r[e][e] = 0));
      for (let e = 0; e < t; e++)
        for (let i = 0; i < t; i++) n[e][i] === 1 && (r[e][i] = 1);
      for (let e = 0; e < t; e++)
        for (let n = 0; n < t; n++)
          for (let i = 0; i < t; i++)
            r[n][i] > r[n][e] + r[e][i] && (r[n][i] = r[n][e] + r[e][i]);
      return r;
    }
    getAdjacencyList() {
      let e = this.vertices.length,
        t = Array(e);
      for (let n = 0; n < e; n++) {
        t[n] = [];
        for (let r = 0; r < e; r++)
          n !== r &&
            this.hasEdge(this.vertices[n].id, this.vertices[r].id) &&
            t[n].push(r);
      }
      return t;
    }
    getSubgraphAdjacencyList(e) {
      let t = e.length,
        n = Array(t);
      for (let r = 0; r < t; r++) {
        n[r] = [];
        for (let i = 0; i < t; i++)
          r !== i && this.hasEdge(e[r], e[i]) && n[r].push(i);
      }
      return n;
    }
    getBridgedRingPerimeter(e, t, n) {
      let r = new Set(t.insiders || []),
        i = e.filter((e) => !r.has(e));
      if (i.length < 3) return i;
      let a = new Set(i),
        o = new Map();
      for (let e = 0; e < i.length; e++) {
        let t = i[e],
          n = this.vertices[t].neighbours.filter((e) => a.has(e));
        if ((o.set(t, n), n.length !== 2)) return i;
      }
      let s = a.has(n) ? n : i[0],
        c = [s],
        l = null,
        u = s;
      for (; c.length < i.length;) {
        let e = o.get(u).filter((e) => e !== l);
        if (
          (l === null &&
            e.sort(
              (e, t) =>
                this.vertices[e].value.smilesOrder -
                this.vertices[t].value.smilesOrder,
            ),
          e.length === 0)
        )
          return i;
        let t = e[0];
        if (t === s) return i;
        (c.push(t), (l = u), (u = t));
      }
      return o.get(u).includes(s) ? c : i;
    }
    getBridges() {
      let e = this.vertices.length,
        t = Array(e),
        n = Array(e),
        r = Array(e),
        i = Array(e),
        a = this.getAdjacencyList(),
        o = [];
      (t.fill(!1), i.fill(null), (this._time = 0));
      for (let s = 0; s < e; s++) t[s] || this._bridgeDfs(s, t, n, r, i, a, o);
      return o;
    }
    traverseBF(e, t) {
      let n = this.vertices.length,
        r = Array(n);
      r.fill(!1);
      let i = [e];
      for (; i.length > 0;) {
        let e = i.shift(),
          n = this.vertices[e];
        t(n);
        for (let e = 0; e < n.neighbours.length; e++) {
          let t = n.neighbours[e];
          r[t] || ((r[t] = !0), i.push(t));
        }
      }
    }
    getTreeDepth(e, t) {
      if (e === null || t === null) return 0;
      let n = this.vertices[e].getSpanningTreeNeighbours(t),
        r = 0;
      for (let t = 0; t < n.length; t++) {
        let i = n[t],
          a = this.getTreeDepth(i, e);
        a > r && (r = a);
      }
      return r + 1;
    }
    traverseTree(e, t, n, r = 999999, i = !1, a = 1, o = null) {
      if (
        (o === null && (o = new Uint8Array(this.vertices.length)),
        a > r + 1 || o[e] === 1)
      )
        return;
      o[e] = 1;
      let s = this.vertices[e],
        c = s.getNeighbours(t);
      (!i || a > 1) && n(s);
      for (let t = 0; t < c.length; t++)
        this.traverseTree(c[t], e, n, r, i, a + 1, o);
    }
    static mdsLayout(t, n, r) {
      let i = Array(n);
      for (let e = 0; e < n; e++) {
        i[e] = new Float64Array(n);
        for (let r = 0; r < n; r++) i[e][r] = t[e][r] * t[e][r];
      }
      let a = new Float64Array(n),
        o = 0;
      for (let e = 0; e < n; e++) {
        let t = 0;
        for (let r = 0; r < n; r++) t += i[e][r];
        ((a[e] = t / n), (o += t));
      }
      o /= n * n;
      let s = Array(n);
      for (let e = 0; e < n; e++) {
        s[e] = new Float64Array(n);
        for (let t = 0; t < n; t++)
          s[e][t] = -0.5 * (i[e][t] - a[e] - a[t] + o);
      }
      let c = e._powerIteration(s, n, 200),
        l = e._rayleigh(s, c, n);
      for (let e = 0; e < n; e++)
        for (let t = 0; t < n; t++) s[e][t] -= l * c[e] * c[t];
      let u = e._powerIteration(s, n, 200),
        d = e._rayleigh(s, u, n),
        f = l > 0 ? Math.sqrt(l) : 0,
        p = d > 0 ? Math.sqrt(d) : 0,
        m = new Float32Array(n),
        h = new Float32Array(n);
      for (let e = 0; e < n; e++) ((m[e] = c[e] * f), (h[e] = u[e] * p));
      let g = 0;
      for (let e = 0; e < n; e++)
        for (let t = e + 1; t < n; t++) {
          let n = m[e] - m[t],
            r = h[e] - h[t],
            i = n * n + r * r;
          i > g && (g = i);
        }
      let _ = Math.sqrt(g);
      if (_ > 0) {
        let e = (r * Math.max(2, Math.sqrt(n))) / _;
        for (let t = 0; t < n; t++) ((m[t] *= e), (h[t] *= e));
      }
      return { xs: m, ys: h };
    }
    static _powerIteration(e, t, n) {
      let r = new Float64Array(t);
      for (let e = 0; e < t; e++) r[e] = 1 + 0.1 * e;
      let i = new Float64Array(t);
      for (let a = 0; a < n; a++) {
        for (let n = 0; n < t; n++) {
          let a = 0;
          for (let i = 0; i < t; i++) a += e[n][i] * r[i];
          i[n] = a;
        }
        let n = 0;
        for (let e = 0; e < t; e++) n += i[e] * i[e];
        if (((n = Math.sqrt(n)), n < 1e-12)) break;
        for (let e = 0; e < t; e++) r[e] = i[e] / n;
      }
      return r;
    }
    static _rayleigh(e, t, n) {
      let r = 0;
      for (let i = 0; i < n; i++) {
        let a = 0;
        for (let r = 0; r < n; r++) a += e[i][r] * t[r];
        r += t[i] * a;
      }
      return r;
    }
    kkLayout(t, n, r, i, a, o = 0.1, s = 0.1, c = 2e3, l = 50, u = 1e9) {
      let d = a,
        f = this.getSubgraphDistanceMatrix(t),
        p = t.length,
        m = new Float32Array(p),
        g = new Float32Array(p),
        v = Array(p),
        y = _.polyCircumradius(a, p),
        b = _.centralAngle(p),
        x = 0,
        S = new Set(i.insiders || []),
        C = t.filter((e) => S.has(e)),
        w = !1;
      for (let e = 0; e < p; e++)
        if (this.vertices[t[e]].positioned) {
          w = !0;
          break;
        }
      if (C.length > 0 && C.length <= 2) {
        let e = this.getBridgedRingPerimeter(t, i, r),
          o = new Set(e),
          s = e.length,
          c = s > 2 ? _.polyCircumradius(a, s) : y,
          l = s > 0 ? (Math.PI * 2) / s : b,
          u = 0,
          d = new Map();
        if (o.has(r)) {
          let t = e.indexOf(r),
            i = this.vertices[r];
          i.positioned && (u = h.subtract(i.position, n).angle() - t * l);
        }
        for (let t = 0; t < e.length; t++) {
          let r = e[t],
            i = this.vertices[r];
          if (i.positioned) d.set(r, i.position.clone());
          else {
            let e = new h(
              n.x + Math.cos(u + t * l) * c,
              n.y + Math.sin(u + t * l) * c,
            );
            d.set(r, e);
          }
        }
        let f = Math.max(a * 0.35, c * 0.35);
        for (let e = 0; e < C.length; e++) {
          let t = C[e],
            r = this.vertices[t];
          if (r.positioned) {
            d.set(t, r.position.clone());
            continue;
          }
          let i = r.neighbours.filter((e) => o.has(e)),
            a = new h(n.x, n.y);
          if (i.length > 0) {
            a = new h(0, 0);
            for (let e = 0; e < i.length; e++) a.add(d.get(i[e]));
            a.divide(i.length);
            let t = h.subtract(a, n);
            (t.lengthSq() < 1e-4
              ? (t = new h(Math.cos(e * l), Math.sin(e * l)))
              : t.normalize(),
              (a = t.multiplyScalar(f).add(n.clone())));
          } else ((a.x += Math.cos(e * l) * f), (a.y += Math.sin(e * l) * f));
          d.set(t, a);
        }
        for (var T = p; T--;) {
          let e = this.vertices[t[T]];
          if (e.positioned) ((m[T] = e.position.x), (g[T] = e.position.y));
          else {
            let t = d.get(e.id);
            t
              ? ((m[T] = t.x), (g[T] = t.y))
              : ((m[T] = n.x + Math.cos(x) * y),
                (g[T] = n.y + Math.sin(x) * y));
          }
          ((v[T] = e.positioned), (x += b));
        }
      } else if (!w && p >= 6) {
        let t = e.mdsLayout(f, p, a);
        for (let e = 0; e < p; e++)
          ((m[e] = n.x + t.xs[e]), (g[e] = n.y + t.ys[e]), (v[e] = !1));
      } else
        for (var T = p; T--;) {
          let e = this.vertices[t[T]];
          (e.positioned
            ? ((m[T] = e.position.x), (g[T] = e.position.y))
            : ((m[T] = n.x + Math.cos(x) * y), (g[T] = n.y + Math.sin(x) * y)),
            (v[T] = e.positioned),
            (x += b));
        }
      let E = Array(p);
      for (T = p; T--;) {
        E[T] = Array(p);
        let e = p;
        for (; e--;) E[T][e] = a * f[T][e];
      }
      let D = Array(p);
      for (T = p; T--;) {
        D[T] = Array(p);
        let e = p;
        for (; e--;) D[T][e] = d * f[T][e] ** -2;
      }
      let O = Array(p),
        k = new Float32Array(p),
        A = new Float32Array(p);
      for (T = p; T--;) O[T] = Array(p);
      for (T = p; T--;) {
        let e = m[T],
          t = g[T],
          n = 0,
          r = 0,
          i = p;
        for (; i--;) {
          if (T === i) continue;
          let a = m[i],
            o = g[i],
            s = 1 / Math.sqrt((e - a) * (e - a) + (t - o) * (t - o));
          ((O[T][i] = [
            D[T][i] * (e - a - E[T][i] * (e - a) * s),
            D[T][i] * (t - o - E[T][i] * (t - o) * s),
          ]),
            (O[i][T] = O[T][i]),
            (n += O[T][i][0]),
            (r += O[T][i][1]));
        }
        ((k[T] = n), (A[T] = r));
      }
      let ee = function (e) {
          return [k[e] * k[e] + A[e] * A[e], k[e], A[e]];
        },
        j = function () {
          let e = 0,
            t = 0,
            n = 0,
            r = 0;
          for (T = p; T--;) {
            let [i, a, o] = ee(T);
            i > e && v[T] === !1 && ((e = i), (t = T), (n = a), (r = o));
          }
          return [t, e, n, r];
        },
        te = function (e, t, n) {
          let r = 0,
            i = 0,
            o = 0,
            s = m[e],
            c = g[e],
            l = E[e],
            u = D[e];
          for (T = p; T--;) {
            if (T === e) continue;
            let t = m[T],
              n = g[T],
              a = l[T],
              d = u[T],
              f = (s - t) * (s - t),
              p = 1 / (f + (c - n) * (c - n)) ** 1.5;
            ((r += d * (1 - a * (c - n) * (c - n) * p)),
              (i += d * (1 - a * f * p)),
              (o += d * (a * (s - t) * (c - n) * p)));
          }
          (r === 0 && (r = 0.1), i === 0 && (i = 0.1), o === 0 && (o = 0.1));
          let d = o / r - i / o,
            f,
            h;
          Math.abs(d) < 1e-6
            ? ((h = -t * 0.1), (f = -n * 0.1))
            : ((f = (t / r + n / o) / d), (h = -(o * f + t) / r));
          let _ = Math.sqrt(h * h + f * f);
          if (_ > a) {
            let e = a / _;
            ((h *= e), (f *= e));
          }
          ((m[e] += h), (g[e] += f));
          let v = O[e];
          for (t = 0, n = 0, s = m[e], c = g[e], T = p; T--;) {
            if (e === T) continue;
            let r = m[T],
              i = g[T],
              a = v[T][0],
              o = v[T][1],
              d = 1 / Math.sqrt((s - r) * (s - r) + (c - i) * (c - i));
            ((h = u[T] * (s - r - l[T] * (s - r) * d)),
              (f = u[T] * (c - i - l[T] * (c - i) * d)),
              (v[T] = [h, f]),
              (t += h),
              (n += f),
              (k[T] += h - a),
              (A[T] += f - o));
          }
          ((k[e] = t), (A[e] = n));
        },
        M = 0,
        N = 0,
        P = 0,
        ne = 0,
        F = 0,
        I = 0;
      for (; u > o && c > F;)
        for (F++, [M, u, N, P] = j(), ne = u, I = 0; ne > s && l > I;)
          (I++, te(M, N, P), ([ne, N, P] = ee(M)));
      for (T = p; T--;) {
        let e = t[T],
          n = this.vertices[e];
        ((n.position.x = m[T]),
          (n.position.y = g[T]),
          (n.positioned = !0),
          (n.forcePositioned = !0));
      }
    }
    _bridgeDfs(e, t, n, r, i, a, o) {
      ((t[e] = !0), (n[e] = r[e] = ++this._time));
      for (let s = 0; s < a[e].length; s++) {
        let c = a[e][s];
        t[c]
          ? c !== i[e] && (r[e] = Math.min(r[e], n[c]))
          : ((i[c] = e),
            this._bridgeDfs(c, t, n, r, i, a, o),
            (r[e] = Math.min(r[e], r[c])),
            r[c] > n[e] && o.push([e, c]));
      }
    }
    static getConnectedComponents(t) {
      let n = t.length,
        r = Array(n),
        i = [];
      r.fill(!1);
      for (let a = 0; a < n; a++)
        if (!r[a]) {
          let n = [];
          ((r[a] = !0),
            n.push(a),
            e._ccGetDfs(a, r, t, n),
            n.length > 1 && i.push(n));
        }
      return i;
    }
    static getConnectedComponentCount(t) {
      let n = t.length,
        r = Array(n),
        i = 0;
      r.fill(!1);
      for (let a = 0; a < n; a++)
        r[a] || ((r[a] = !0), i++, e._ccCountDfs(a, r, t));
      return i;
    }
    static _ccCountDfs(t, n, r) {
      for (let i = 0; i < r[t].length; i++)
        !r[t][i] || n[i] || t === i || ((n[i] = !0), e._ccCountDfs(i, n, r));
    }
    static _ccGetDfs(t, n, r, i) {
      for (let a = 0; a < r[t].length; a++)
        !r[t][a] ||
          n[a] ||
          t === a ||
          ((n[a] = !0), i.push(a), e._ccGetDfs(a, n, r, i));
    }
  },
  O = class e {
    static extend() {
      let t = {},
        n = !1,
        r = 0,
        i = arguments.length;
      Object.prototype.toString.call(arguments[0]) === `[object Boolean]` &&
        ((n = arguments[0]), r++);
      let a = function (r) {
        for (let i in r)
          Object.prototype.hasOwnProperty.call(r, i) &&
            (n && Object.prototype.toString.call(r[i]) === `[object Object]`
              ? (t[i] = e.extend(!0, t[i], r[i]))
              : (t[i] = r[i]));
      };
      for (; r < i; r++) {
        let e = arguments[r];
        a(e);
      }
      return t;
    }
  },
  k = class e {
    constructor(e) {
      this.size = e;
      let t = e === 0 ? 0 : ((e - 1) >>> 5) + 1;
      this.words = new Uint32Array(t);
    }
    set(e) {
      this.words[e >>> 5] |= 1 << (e & 31);
    }
    get(e) {
      return !!(this.words[e >>> 5] & (1 << (e & 31)));
    }
    _assertSameSize(e, t) {
      if (this.size !== e.size)
        throw Error(
          `BitSet.` +
            t +
            `: size mismatch (` +
            this.size +
            ` vs ` +
            e.size +
            `). Operands must share the same bit length.`,
        );
    }
    xor(t) {
      this._assertSameSize(t, `xor`);
      let n = new e(this.size);
      for (let e = 0; e < this.words.length; e++)
        n.words[e] = this.words[e] ^ t.words[e];
      return n;
    }
    and(t) {
      this._assertSameSize(t, `and`);
      let n = new e(this.size);
      for (let e = 0; e < this.words.length; e++)
        n.words[e] = this.words[e] & t.words[e];
      return n;
    }
    or(t) {
      this._assertSameSize(t, `or`);
      let n = new e(this.size);
      for (let e = 0; e < this.words.length; e++)
        n.words[e] = this.words[e] | t.words[e];
      return n;
    }
    orWith(e) {
      this._assertSameSize(e, `orWith`);
      for (let t = 0; t < this.words.length; t++) this.words[t] |= e.words[t];
    }
    isSubsetOf(e) {
      this._assertSameSize(e, `isSubsetOf`);
      for (let t = 0; t < this.words.length; t++)
        if ((this.words[t] & ~e.words[t]) !== 0) return !1;
      return !0;
    }
    clone() {
      let t = new e(this.size);
      for (let e = 0; e < this.words.length; e++) t.words[e] = this.words[e];
      return t;
    }
    isEmpty() {
      for (let e = 0; e < this.words.length; e++)
        if (this.words[e] !== 0) return !1;
      return !0;
    }
    popcount() {
      let e = 0;
      for (let t = 0; t < this.words.length; t++) {
        let n = this.words[t];
        for (; n !== 0;) ((n &= n - 1), e++);
      }
      return e;
    }
    clz() {
      for (let e = this.words.length - 1; e >= 0; e--)
        if (this.words[e] !== 0) {
          let t = e * 32 + (31 - Math.clz32(this.words[e]));
          return this.size - 1 - t;
        }
      return this.size;
    }
  },
  A = 2147483647;
function ee(e, t, n, r) {
  let i = e.length,
    a = new Int32Array(i).fill(A),
    o = Array(i).fill(!1),
    s = Array(i).fill(null),
    c = new Int32Array(i);
  ((a[t] = 0), (o[t] = !0), (c[t] = 1), (s[t] = { type: `source`, vertex: t }));
  let l = [t],
    u = 1;
  for (let i = 0; i < u; i++) {
    let d = l[i],
      f = a[d] + 1;
    if (!(f > n))
      for (let n = 0; n < e[d].length; n++) {
        let i = e[d][n];
        if (f < a[i])
          ((a[i] = f),
            (s[i] = { type: `seq`, parent: s[d], vertex: i, dist: f }),
            (o[i] = o[d] && r[i] < r[t]),
            (c[i] = c[d]),
            (l[u++] = i));
        else if (f === a[i]) {
          if (!(o[d] && r[i] < r[t])) continue;
          let e = { type: `seq`, parent: s[d], vertex: i, dist: f };
          o[i]
            ? ((s[i] = { type: `branch`, left: s[i], right: e }),
              (c[i] += c[d]))
            : ((o[i] = !0), (s[i] = e), (c[i] = c[d]));
        }
      }
  }
  function d(e, t) {
    if (!e) return [];
    let n = Array(t),
      r = e;
    for (; r;)
      if (r.type === `source`) {
        n[0] = r.vertex;
        break;
      } else
        r.type === `seq`
          ? ((n[r.dist] = r.vertex), (r = r.parent))
          : (r = r.left);
    return n;
  }
  function f(e, t) {
    if (!e) return [];
    if (e.type === `source`) {
      let n = Array(t);
      return ((n[0] = e.vertex), [n]);
    }
    if (e.type === `seq`) {
      let n = f(e.parent, t);
      for (let t = 0; t < n.length; t++) n[t][e.dist] = e.vertex;
      return n;
    }
    {
      let n = f(e.left, t),
        r = f(e.right, t);
      return n.concat(r);
    }
  }
  return {
    distTo: a,
    nPathsTo: c,
    precedes: o,
    pathTo(e) {
      return e < 0 || e >= i || !s[e] ? [] : d(s[e], a[e] + 1);
    },
    pathsTo(e) {
      return e < 0 || e >= i || !s[e] ? [] : f(s[e], a[e] + 1);
    },
    isPrecedingPathTo(e) {
      return e >= 0 && e < i && o[e];
    },
  };
}
function j(e, t) {
  let n = e.length;
  for (let r = 1; r < n; r++) if (e[r] === t[r]) return !1;
  return !0;
}
function te(e, t, n) {
  let r = new k(n),
    i = e.length - 1;
  for (let n = 0; n < i; n++) {
    let i = e[n],
      a = e[n + 1],
      o = i < a ? i + `,` + a : a + `,` + i;
    r.set(t.get(o));
  }
  return r;
}
function M(e) {
  let t = new Map(),
    n = e.length;
  for (let r = 0; r < n; r++)
    for (let n = 0; n < e[r].length; n++) {
      let i = e[r][n];
      i > r && t.set(r + `,` + i, t.size);
    }
  return t;
}
function N(e) {
  let t = e.length,
    n = 0,
    r = Array(t).fill(0);
  for (let r = 0; r < t; r++) e[r].length > n && (n = e[r].length);
  let i = Array(n + 1);
  for (let e = 0; e <= n; e++) i[e] = [];
  for (let n = 0; n < t; n++) i[e[n].length].push(n);
  let a = 0;
  for (let e = 0; e <= n; e++) for (let t of i[e]) r[t] = a++;
  return r;
}
function P(e) {
  let t = e.length,
    n = N(e),
    r = new Map(),
    i = M(e),
    a = i.size,
    o = Array(t);
  for (let e = 0; e < t; e++) o[n[e]] = e;
  function s(e) {
    let t = te(e, i, a),
      n = e.length - 1,
      o = r.get(n);
    (o || ((o = []), r.set(n, o)), o.push({ path: e, edgeVector: t }));
  }
  let c = Array(t);
  for (let r = 2; r < t; r++) {
    let i = o[r],
      a = ee(e, i, Math.floor(t / 2), n);
    for (let t = 0; t < r; t++) {
      let r = o[t];
      if (!a.isPrecedingPathTo(r)) continue;
      let i = 0;
      for (let t = 0; t < e[r].length; t++) {
        let o = e[r][t];
        if (!a.isPrecedingPathTo(o)) continue;
        let l = a.distTo[o],
          u = a.distTo[r];
        if (l + 1 === u) c[i++] = o;
        else if (l === u && n[o] < n[r]) {
          let e = a.pathTo(r),
            t = a.pathTo(o);
          j(t, e) && s(e.concat(t.reverse()));
        }
      }
      for (let e = 0; e < i; e++)
        for (let t = e + 1; t < i; t++) {
          let n = a.pathTo(c[e]),
            i = a.pathTo(c[t]);
          j(n, i) && s(n.concat([r], i.reverse()));
        }
    }
  }
  return { cycles: r, edgeIndex: i, nEdges: a };
}
var ne = class {
    constructor(e, t) {
      ((this._n = e),
        (this._max = t),
        (this._rows = Array(t)),
        (this._indices = new Int32Array(t)),
        (this._m = 0));
    }
    add(e) {
      ((this._rows[this._m] = e),
        (this._indices[this._m] = this._m),
        this._m++);
    }
    swap(e, t) {
      let n = this._rows[e],
        r = this._indices[e];
      ((this._rows[e] = this._rows[t]),
        (this._indices[e] = this._indices[t]),
        (this._rows[t] = n),
        (this._indices[t] = r));
    }
    _rowIndex(e) {
      for (let t = 0; t < this._m; t++) if (this._indices[t] === e) return t;
      return -1;
    }
    eliminated(e) {
      return this._rows[this._rowIndex(e)].isEmpty();
    }
    eliminate() {
      return this._eliminate(0, 0);
    }
    _eliminate(e, t) {
      for (; e < this._n && t < this._m;) {
        let n = -1;
        for (let r = t; r < this._m; r++)
          if (this._rows[r].get(e)) {
            n = r;
            break;
          }
        if (n < 0) return this._eliminate(e + 1, t);
        n !== t && this.swap(n, t);
        for (let n = t + 1; n < this._m; n++)
          this._rows[n].get(e) &&
            (this._rows[n] = this._rows[n].xor(this._rows[t]));
        t++;
      }
      return t;
    }
  },
  F = class {
    constructor(e) {
      ((this._members = []),
        (this._edgesOfBasis = new k(e)),
        (this._nEdges = e));
    }
    add(e) {
      (this._members.push(e), this._edgesOfBasis.orWith(e.edgeVector));
    }
    members() {
      return this._members;
    }
    size() {
      return this._members.length;
    }
    isSubsetOfBasis(e) {
      return e.edgeVector.isSubsetOf(this._edgesOfBasis);
    }
    isIndependent(e) {
      if (this._members.length === 0 || !this.isSubsetOfBasis(e)) return !0;
      let t = new ne(this._nEdges, this._members.length + 1);
      for (let e = 0; e < this._members.length; e++)
        t.add(this._members[e].edgeVector.clone());
      return (
        t.add(e.edgeVector.clone()),
        t.eliminate(),
        !t.eliminated(this._members.length)
      );
    }
  };
function I(e) {
  let t = e.length;
  if (t < 3) return [];
  let { cycles: n, nEdges: r } = P(e),
    i = new F(r),
    a = r - t + 1,
    o = [...n.keys()].sort((e, t) => e - t);
  for (let e = 0; e < o.length; e++) {
    let t = n.get(o[e]);
    for (let e = 0; e < t.length && !(i.size() >= a); e++)
      i.isIndependent(t[e]) && i.add(t[e]);
  }
  return i.members().map((e) => e.path);
}
function re(e) {
  if (e.length < 3) return [];
  let { cycles: t, nEdges: n } = P(e),
    r = new F(n),
    i = [],
    a = [...t.keys()].sort((e, t) => e - t);
  for (let e = 0; e < a.length; e++) {
    let n = t.get(a[e]),
      o = [];
    for (let e = 0; e < n.length; e++)
      r.isIndependent(n[e]) && (o.push(n[e]), i.push(n[e].path));
    for (let e = 0; e < o.length; e++) r.add(o[e]);
  }
  return i;
}
function ie(e) {
  let t = e.length,
    n = Array(t);
  for (let r = 0; r < t; r++) {
    n[r] = [];
    for (let i = 0; i < t; i++) e[r][i] === 1 && n[r].push(i);
  }
  return n;
}
var ae = class e {
    static getRings(t, n = !1) {
      return e._findRings(t, !1);
    }
    static getRingsForLayout(t, n = !1) {
      return e._findRings(t, !0);
    }
    static _findRings(e, t) {
      let n = e.getComponentsAdjacencyMatrix();
      if (n.length === 0) return [];
      let r = D.getConnectedComponents(n),
        i = [];
      for (let n = 0; n < r.length; n++) {
        let a = r[n],
          o = ie(e.getSubgraphAdjacencyMatrix([...a])),
          s = t ? re(o) : I(o);
        for (let e = 0; e < s.length; e++) {
          let t = s[e],
            n = [];
          for (let e = 0; e < t.length - 1; e++) n.push(a[t[e]]);
          i.push(n);
        }
      }
      return i;
    }
    static matrixToString(e) {
      let t = ``;
      for (let n = 0; n < e.length; n++) {
        for (let r = 0; r < e[n].length; r++) t += e[n][r] + ` `;
        t += `
`;
      }
      return t;
    }
    static areSetsEqual(e, t) {
      if (e.size !== t.size) return !1;
      for (let n of e) if (!t.has(n)) return !1;
      return !0;
    }
    static isSupersetOf(e, t) {
      for (let n of t) if (!e.has(n)) return !1;
      return !0;
    }
  },
  oe = class e {
    constructor(e) {
      ((this.graph = null),
        (this.doubleBondConfigCount = 0),
        (this.doubleBondConfig = null),
        (this.ringIdCounter = 0),
        (this.ringConnectionIdCounter = 0),
        (this.canvasWrapper = null),
        (this.totalOverlapScore = 0),
        (this.defaultOptions = {
          width: 500,
          height: 500,
          scale: 0,
          bondThickness: 1,
          bondLength: 30,
          shortBondLength: 0.8,
          bondSpacing: 0.17 * 30,
          atomVisualization: `default`,
          isomeric: !0,
          debug: !1,
          terminalCarbons: !1,
          showCarbons: `default`,
          explicitHydrogens: !0,
          overlapSensitivity: 0.42,
          overlapResolutionIterations: 1,
          compactDrawing: !0,
          fontFamily: `Arial, Helvetica, sans-serif`,
          fontSizeLarge: 11,
          fontSizeSmall: 3,
          padding: 10,
          experimentalSSSR: !1,
          experimentalWeights: !1,
          kkThreshold: 0.1,
          kkInnerThreshold: 0.1,
          kkMaxIteration: 2e4,
          kkMaxInnerIteration: 50,
          kkMaxEnergy: 1e9,
          weights: {
            colormap: null,
            additionalPadding: 20,
            sigma: 10,
            interval: 0,
            opacity: 1,
          },
          themes: {
            dark: {
              FOREGROUND: `#ffffff`,
              BACKGROUND: `#141414`,
              C: `#ffffff`,
              O: `#e74c3c`,
              N: `#3498db`,
              F: `#27ae60`,
              CL: `#16a085`,
              BR: `#d35400`,
              I: `#8e44ad`,
              P: `#d35400`,
              S: `#f1c40f`,
              B: `#e67e22`,
              SI: `#e67e22`,
              H: `#aaaaaa`,
            },
            light: {
              FOREGROUND: `#222222`,
              BACKGROUND: `#ffffff`,
              C: `#222222`,
              O: `#e74c3c`,
              N: `#3498db`,
              F: `#27ae60`,
              CL: `#16a085`,
              BR: `#d35400`,
              I: `#8e44ad`,
              P: `#d35400`,
              S: `#f1c40f`,
              B: `#e67e22`,
              SI: `#e67e22`,
              H: `#666666`,
            },
            oldschool: {
              FOREGROUND: `#000000`,
              BACKGROUND: `#ffffff`,
              C: `#000000`,
              O: `#000000`,
              N: `#000000`,
              F: `#000000`,
              CL: `#000000`,
              BR: `#000000`,
              I: `#000000`,
              P: `#000000`,
              S: `#000000`,
              B: `#000000`,
              SI: `#000000`,
              H: `#000000`,
            },
            solarized: {
              FOREGROUND: `#586e75`,
              BACKGROUND: `#eee8d5`,
              C: `#586e75`,
              O: `#dc322f`,
              N: `#268bd2`,
              F: `#859900`,
              CL: `#16a085`,
              BR: `#cb4b16`,
              I: `#6c71c4`,
              P: `#d33682`,
              S: `#b58900`,
              B: `#2aa198`,
              SI: `#2aa198`,
              H: `#657b83`,
            },
            "solarized-dark": {
              FOREGROUND: `#93a1a1`,
              BACKGROUND: `#073642`,
              C: `#93a1a1`,
              O: `#dc322f`,
              N: `#268bd2`,
              F: `#859900`,
              CL: `#16a085`,
              BR: `#cb4b16`,
              I: `#6c71c4`,
              P: `#d33682`,
              S: `#b58900`,
              B: `#2aa198`,
              SI: `#2aa198`,
              H: `#839496`,
            },
            matrix: {
              FOREGROUND: `#678c61`,
              BACKGROUND: `#ffffff`,
              C: `#678c61`,
              O: `#2fc079`,
              N: `#4f7e7e`,
              F: `#90d762`,
              CL: `#82d967`,
              BR: `#23755a`,
              I: `#409931`,
              P: `#c1ff8a`,
              S: `#faff00`,
              B: `#50b45a`,
              SI: `#409931`,
              H: `#426644`,
            },
            github: {
              FOREGROUND: `#24292f`,
              BACKGROUND: `#ffffff`,
              C: `#24292f`,
              O: `#cf222e`,
              N: `#0969da`,
              F: `#2da44e`,
              CL: `#6fdd8b`,
              BR: `#bc4c00`,
              I: `#8250df`,
              P: `#bf3989`,
              S: `#d4a72c`,
              B: `#fb8f44`,
              SI: `#bc4c00`,
              H: `#57606a`,
            },
            carbon: {
              FOREGROUND: `#161616`,
              BACKGROUND: `#ffffff`,
              C: `#161616`,
              O: `#da1e28`,
              N: `#0f62fe`,
              F: `#198038`,
              CL: `#007d79`,
              BR: `#fa4d56`,
              I: `#8a3ffc`,
              P: `#ff832b`,
              S: `#f1c21b`,
              B: `#8a3800`,
              SI: `#e67e22`,
              H: `#525252`,
            },
            cyberpunk: {
              FOREGROUND: `#ea00d9`,
              BACKGROUND: `#ffffff`,
              C: `#ea00d9`,
              O: `#ff3131`,
              N: `#0abdc6`,
              F: `#00ff9f`,
              CL: `#00fe00`,
              BR: `#fe9f20`,
              I: `#ff00ff`,
              P: `#fe7f00`,
              S: `#fcee0c`,
              B: `#ff00ff`,
              SI: `#ffffff`,
              H: `#913cb1`,
            },
            gruvbox: {
              FOREGROUND: `#665c54`,
              BACKGROUND: `#fbf1c7`,
              C: `#665c54`,
              O: `#cc241d`,
              N: `#458588`,
              F: `#98971a`,
              CL: `#79740e`,
              BR: `#d65d0e`,
              I: `#b16286`,
              P: `#af3a03`,
              S: `#d79921`,
              B: `#689d6a`,
              SI: `#427b58`,
              H: `#7c6f64`,
            },
            "gruvbox-dark": {
              FOREGROUND: `#ebdbb2`,
              BACKGROUND: `#282828`,
              C: `#ebdbb2`,
              O: `#cc241d`,
              N: `#458588`,
              F: `#98971a`,
              CL: `#b8bb26`,
              BR: `#d65d0e`,
              I: `#b16286`,
              P: `#fe8019`,
              S: `#d79921`,
              B: `#8ec07c`,
              SI: `#83a598`,
              H: `#bdae93`,
            },
            custom: {
              FOREGROUND: `#222222`,
              BACKGROUND: `#ffffff`,
              C: `#222222`,
              O: `#e74c3c`,
              N: `#3498db`,
              F: `#27ae60`,
              CL: `#16a085`,
              BR: `#d35400`,
              I: `#8e44ad`,
              P: `#d35400`,
              S: `#f1c40f`,
              B: `#e67e22`,
              SI: `#e67e22`,
              H: `#666666`,
            },
          },
        }),
        (this.opts = O.extend(!0, this.defaultOptions, e)),
        [`none`, `default`, `terminal`, `acyclic`, `all`].indexOf(
          this.opts.showCarbons,
        ) === -1 && (this.opts.showCarbons = `default`),
        (this.opts.halfBondSpacing = this.opts.bondSpacing / 2),
        (this.opts.bondLengthSq = this.opts.bondLength * this.opts.bondLength),
        (this.opts.halfFontSizeLarge = this.opts.fontSizeLarge / 2),
        (this.opts.quarterFontSizeLarge = this.opts.fontSizeLarge / 4),
        (this.opts.fifthFontSizeSmall = this.opts.fontSizeSmall / 5),
        (this.theme = this.opts.themes.dark));
    }
    static getEffectiveShowCarbonsMode(e) {
      let t = [`none`, `default`, `terminal`, `acyclic`, `all`],
        n = e.showCarbons;
      return (
        (n == null || t.indexOf(n) === -1) && (n = `default`),
        n === "default" && e.terminalCarbons ? `terminal` : n
      );
    }
    draw(e, t, n = `light`, r = !1) {
      (this.initDraw(e, n, r),
        this.infoOnly ||
          ((this.themeManager = new x(this.opts.themes, n)),
          (this.canvasWrapper = new C(t, this.themeManager, this.opts))),
        r ||
          (this.processGraph(),
          this.canvasWrapper.scale(this.graph.vertices),
          this.drawEdges(this.opts.debug),
          this.drawVertices(this.opts.debug),
          this.canvasWrapper.reset(),
          this.opts.debug &&
            console.debug(`DrawerBase::draw()`, {
              graph: this.graph,
              rings: this.rings,
              ringConnections: this.ringConnections,
            })));
    }
    edgeRingCount(e) {
      let t = this.graph.edges[e],
        n = this.graph.vertices[t.sourceId],
        r = this.graph.vertices[t.targetId];
      return Math.min(n.value.rings.length, r.value.rings.length);
    }
    getBridgedRings() {
      return this.rings.filter((e) => e.isBridged);
    }
    getFusedRings() {
      return this.rings.filter((e) => e.isFused);
    }
    getSpiros() {
      return this.rings.filter((e) => e.isSpiro);
    }
    printRingInfo() {
      let e = ``;
      for (let t = 0; t < this.rings.length; t++) {
        let n = this.rings[t];
        ((e += n.id + `;`),
          (e += n.members.length + `;`),
          (e += n.neighbours.length + `;`),
          (e += n.isSpiro ? `true;` : `false;`),
          (e += n.isFused ? `true;` : `false;`),
          (e += n.isBridged ? `true;` : `false;`),
          (e += n.rings.length + `;`),
          (e += `
`));
      }
      return e;
    }
    rotateDrawing() {
      let e = 0,
        t = 0,
        n = 0;
      for (let r = 0; r < this.graph.vertices.length; r++) {
        let i = this.graph.vertices[r];
        if (i.value.isDrawn)
          for (let a = r + 1; a < this.graph.vertices.length; a++) {
            let o = this.graph.vertices[a];
            if (!o.value.isDrawn) continue;
            let s = i.position.distanceSq(o.position);
            s > n && ((n = s), (e = r), (t = a));
          }
      }
      let r = -h
        .subtract(
          this.graph.vertices[e].position,
          this.graph.vertices[t].position,
        )
        .angle();
      if (!isNaN(r)) {
        let e = r % 0.523599;
        e < 0.2617995 ? (r -= e) : (r += 0.523599 - e);
        for (let e = 0; e < this.graph.vertices.length; e++)
          e !== t &&
            this.graph.vertices[e].position.rotateAround(
              r,
              this.graph.vertices[t].position,
            );
        for (let e = 0; e < this.rings.length; e++)
          this.rings[e].center.rotateAround(r, this.graph.vertices[t].position);
      }
    }
    getTotalOverlapScore() {
      return this.totalOverlapScore;
    }
    getRingCount() {
      return this.rings.length;
    }
    hasBridgedRing() {
      return this.bridgedRing;
    }
    getHeavyAtomCount() {
      let e = 0;
      for (let t = 0; t < this.graph.vertices.length; t++)
        this.graph.vertices[t].value.element !== `H` && e++;
      return e;
    }
    getMolecularFormula(e = null) {
      let t = ``,
        n = new Map(),
        r = e === null ? this.graph : new D(e, this.opts.isomeric);
      for (let e = 0; e < r.vertices.length; e++) {
        let t = r.vertices[e].value,
          i = n.get(t.element) || 0;
        n.set(t.element, i + 1);
        let a = t.countImplicitHydrogens();
        if (a) {
          let e = n.get(`H`) || 0;
          n.set(`H`, e + a);
        }
      }
      if (n.has(`C`)) {
        let e = n.get(`C`);
        ((t += `C` + (e > 1 ? e : ``)), n.delete(`C`));
      }
      if (n.has(`H`)) {
        let e = n.get(`H`);
        ((t += `H` + (e > 1 ? e : ``)), n.delete(`H`));
      }
      return (
        Object.keys(m.atomicNumbers)
          .sort()
          .map((e) => {
            if (n.has(e)) {
              let r = n.get(e);
              t += e + (r > 1 ? r : ``);
            }
          }),
        t
      );
    }
    getRingbondType(e, t) {
      if (e.value.getRingbondCount() < 1 || t.value.getRingbondCount() < 1)
        return null;
      for (let n = 0; n < e.value.ringbonds.length; n++)
        for (let r = 0; r < t.value.ringbonds.length; r++)
          if (e.value.ringbonds[n].id === t.value.ringbonds[r].id)
            return e.value.ringbonds[n].bondType === `-`
              ? t.value.ringbonds[r].bond
              : e.value.ringbonds[n].bond;
      return null;
    }
    initDraw(e, t, n, r) {
      ((this.data = e),
        (this.infoOnly = n),
        (this.ringIdCounter = 0),
        (this.ringConnectionIdCounter = 0),
        (this.graph = new D(e, this.opts.isomeric)),
        (this.rings = []),
        (this.ringConnections = []),
        (this.originalRings = []),
        (this.originalRingConnections = []),
        (this.bridgedRing = !1),
        (this.doubleBondConfigCount = null),
        (this.doubleBondConfig = null),
        (this.highlight_atoms = r),
        this.initRings(),
        this.initHydrogens());
    }
    processGraph() {
      (this.position(),
        this.fixDoubleBondStereo(),
        this.restoreRingInformation(),
        this.resolvePrimaryOverlaps());
      let e = this.getOverlapScore();
      this.totalOverlapScore = this.getOverlapScore().total;
      for (let t = 0; t < this.opts.overlapResolutionIterations; t++)
        for (let t = 0; t < this.graph.edges.length; t++) {
          let n = this.graph.edges[t];
          if (this.isEdgeRotatable(n)) {
            let t = this.graph.getTreeDepth(n.sourceId, n.targetId),
              r = this.graph.getTreeDepth(n.targetId, n.sourceId),
              i = n.targetId,
              a = n.sourceId;
            if (
              (t > r && ((i = n.sourceId), (a = n.targetId)),
              this.getSubtreeOverlapScore(a, i, e.vertexScores).value >
                this.opts.overlapSensitivity)
            ) {
              let t = this.graph.vertices[i],
                n = this.graph.vertices[a],
                r = n.getNeighbours(i);
              if (r.length === 1) {
                let e = this.graph.vertices[r[0]],
                  i = e.position.getRotateAwayFromAngle(
                    t.position,
                    n.position,
                    _.toRad(120),
                  );
                this.rotateSubtree(e.id, n.id, i, n.position);
                let a = this.getOverlapScore().total;
                a > this.totalOverlapScore
                  ? this.rotateSubtree(e.id, n.id, -i, n.position)
                  : (this.totalOverlapScore = a);
              } else if (r.length === 2) {
                if (n.value.rings.length !== 0 && t.value.rings.length !== 0)
                  continue;
                let e = this.graph.vertices[r[0]],
                  i = this.graph.vertices[r[1]];
                if (e.value.rings.length === 1 && i.value.rings.length === 1) {
                  if (e.value.rings[0] !== i.value.rings[0]) continue;
                  let r = e.value.rings[0];
                  if (this.getRingExternalConnectionCount(r) !== 1) continue;
                  let a = 0,
                    o = this.totalOverlapScore,
                    s = this.getRing(r),
                    c = _.centralAngle(s.getSize()),
                    l = Math.max(1, Math.floor(s.getSize() / 2));
                  for (let e = 1; e <= l; e++) {
                    let r = c * e;
                    this.rotateSubtree(n.id, t.id, r, n.position);
                    let i = this.getOverlapScore().total;
                    (i < o && ((o = i), (a = r)),
                      this.rotateSubtree(n.id, t.id, -r, n.position),
                      this.rotateSubtree(n.id, t.id, -r, n.position),
                      (i = this.getOverlapScore().total),
                      i < o && ((o = i), (a = -r)),
                      this.rotateSubtree(n.id, t.id, r, n.position));
                  }
                  a !== 0 &&
                    (this.rotateSubtree(n.id, t.id, a, n.position),
                    (this.totalOverlapScore = o));
                } else {
                  if (e.value.rings.length !== 0 || i.value.rings.length !== 0)
                    continue;
                  {
                    let r = e.position.getRotateAwayFromAngle(
                        t.position,
                        n.position,
                        _.toRad(120),
                      ),
                      a = i.position.getRotateAwayFromAngle(
                        t.position,
                        n.position,
                        _.toRad(120),
                      );
                    (this.rotateSubtree(e.id, n.id, r, n.position),
                      this.rotateSubtree(i.id, n.id, a, n.position));
                    let o = this.getOverlapScore().total;
                    o > this.totalOverlapScore
                      ? (this.rotateSubtree(e.id, n.id, -r, n.position),
                        this.rotateSubtree(i.id, n.id, -a, n.position))
                      : (this.totalOverlapScore = o);
                  }
                }
              }
              e = this.getOverlapScore();
            }
          }
        }
      (this.resolveSecondaryOverlaps(e.scores),
        this.resolveRigidRingOverlaps(),
        (e = this.getOverlapScore()),
        this.resolveSecondaryOverlaps(e.scores),
        this.opts.isomeric && this.annotateStereochemistry(),
        this.opts.compactDrawing &&
          this.opts.atomVisualization === "default" &&
          this.initPseudoElements(),
        this.rotateDrawing());
    }
    static flipEZ(e) {
      return e === `/` ? `\\` : e === `\\` ? `/` : e;
    }
    static getRingbondType(t, n) {
      return t && t !== `-` ? t : n && n !== `-` ? e.flipEZ(n) : `-`;
    }
    initRings() {
      let t = new Map();
      for (let n = this.graph.vertices.length - 1; n >= 0; n--) {
        let r = this.graph.vertices[n];
        if (r.value.ringbonds.length !== 0)
          for (let n = 0; n < r.value.ringbonds.length; n++) {
            let i = r.value.ringbonds[n].id,
              a = r.value.ringbonds[n].bond;
            if (!t.has(i)) t.set(i, [r.id, a]);
            else {
              let o = r.id,
                s = t.get(i)[0],
                c = t.get(i)[1],
                l = new E(o, s, 1);
              l.setBondType(e.getRingbondType(a, c));
              let u = this.graph.addEdge(l),
                d = this.graph.vertices[s];
              (r.addRingbondChild(s, n),
                r.value.addNeighbouringElement(d.value.element));
              let f = 0;
              for (let e = 0; e < d.value.ringbonds.length; e++)
                if (d.value.ringbonds[e].id === i) {
                  f = e;
                  break;
                }
              (d.addRingbondChild(o, f),
                d.value.addNeighbouringElement(r.value.element),
                r.edges.push(u),
                d.edges.push(u),
                t.delete(i));
            }
          }
      }
      let n = ae.getRings(this.graph, this.opts.experimentalSSSR);
      if (n !== null && n.length !== 0) {
        for (let e = 0; e < n.length; e++) {
          let t = [...n[e]],
            r = this.addRing(new b(t));
          for (let e = 0; e < t.length; e++)
            this.graph.vertices[t[e]].value.rings.push(r);
        }
        for (let e = 0; e < this.rings.length - 1; e++)
          for (let t = e + 1; t < this.rings.length; t++) {
            let n = this.rings[e],
              r = this.rings[t],
              i = new y(n, r);
            i.vertices.size > 0 && this.addRingConnection(i);
          }
        for (let e = 0; e < this.rings.length; e++) {
          let t = this.rings[e];
          t.neighbours = y.getNeighbours(this.ringConnections, t.id);
        }
        for (let e = 0; e < this.rings.length; e++) {
          let t = this.rings[e];
          this.graph.vertices[t.members[0]].value.addAnchoredRing(t.id);
        }
        for (
          this.markCageRingSystems(), this.backupRingInformation();
          this.rings.length > 0;
        ) {
          let e = -1;
          for (let t = 0; t < this.rings.length; t++) {
            let n = this.rings[t];
            this.isPartOfBridgedRing(n.id) && !n.isBridged && (e = n.id);
          }
          if (e === -1) break;
          let t = this.getRing(e),
            n = this.getBridgedRingRings(t.id);
          ((this.bridgedRing = !0),
            this.createBridgedRing(n, t.members[0]),
            (this.bridgedRing = !1));
          for (let e = 0; e < n.length; e++) this.removeRing(n[e]);
        }
      }
    }
    initHydrogens() {
      if (!this.opts.explicitHydrogens)
        for (let e of this.graph.vertices) {
          if (e.value.element !== `H` || e.neighbours.length !== 1) continue;
          let t = this.graph.vertices[e.neighbours[0]];
          (!t.value.isStereoCenter ||
            (t.value.rings.length < 2 && !t.value.bridgedRing) ||
            (t.value.bridgedRing && t.value.originalRings.length < 2)) &&
            (e.value.isDrawn = !1);
        }
    }
    getBridgedRingRings(e) {
      let t = [],
        n = (e) => {
          let r = this.getRing(e);
          t.push(e);
          for (let i = 0; i < r.neighbours.length; i++) {
            let a = r.neighbours[i];
            t.indexOf(a) === -1 &&
              a !== e &&
              y.isBridge(this.ringConnections, this.graph.vertices, e, a) &&
              n(a);
          }
        };
      n(e);
      let r = !0;
      for (; r;) {
        r = !1;
        for (let e = 0; e < this.ringConnections.length; e++) {
          let n = this.ringConnections[e];
          if (n.vertices.size < 2) continue;
          let i = t.indexOf(n.firstRingId) !== -1,
            a = t.indexOf(n.secondRingId) !== -1;
          i && !a
            ? (t.push(n.secondRingId), (r = !0))
            : a && !i && (t.push(n.firstRingId), (r = !0));
        }
      }
      return f.unique(t);
    }
    isPartOfBridgedRing(e) {
      for (let t = 0; t < this.ringConnections.length; t++)
        if (
          this.ringConnections[t].containsRing(e) &&
          this.ringConnections[t].isBridge(this.graph.vertices)
        )
          return !0;
      return !1;
    }
    markCageRingSystems() {
      let e = new Map(),
        t = new Map();
      for (let n = 0; n < this.rings.length; n++)
        (e.set(this.rings[n].id, 0), t.set(this.rings[n].id, []));
      for (let n = 0; n < this.ringConnections.length; n++) {
        let r = this.ringConnections[n];
        r.vertices.size < 2 ||
          (e.set(r.firstRingId, e.get(r.firstRingId) + 1),
          e.set(r.secondRingId, e.get(r.secondRingId) + 1),
          t.get(r.firstRingId).push(r.secondRingId),
          t.get(r.secondRingId).push(r.firstRingId));
      }
      let n = new Set(),
        r = new Set();
      for (let i = 0; i < this.rings.length; i++) {
        let a = this.rings[i].id;
        if (n.has(a)) continue;
        let o = [],
          s = [a];
        for (n.add(a); s.length > 0;) {
          let e = s.shift();
          o.push(e);
          let r = t.get(e);
          for (let e = 0; e < r.length; e++) {
            let t = r[e];
            n.has(t) || (n.add(t), s.push(t));
          }
        }
        if (this.isCageRingComponent(o, e))
          for (let e = 0; e < o.length; e++) r.add(o[e]);
      }
      if (r.size !== 0)
        for (let e = 0; e < this.ringConnections.length; e++) {
          let t = this.ringConnections[e];
          t.vertices.size < 2 ||
            (r.has(t.firstRingId) &&
              r.has(t.secondRingId) &&
              (t.isForcedBridge = !0));
        }
    }
    isCageRingComponent(e, t) {
      let n = 0;
      for (let r = 0; r < e.length; r++) t.get(e[r]) >= 3 && n++;
      if (n < 2) return !1;
      let r = this.getRingSystemStats(e);
      return (
        r.atomCount > 0 &&
        r.edgeCount > 0 &&
        r.nonCageAtomCount === 0 &&
        r.boundaryEdgeRatio <= 0.5
      );
    }
    getRingSystemStats(e) {
      let t = new Map(),
        n = new Set(),
        r = new Map();
      for (let i = 0; i < e.length; i++) {
        let a = this.getRing(e[i]),
          o = new Set(a.members);
        t.set(e[i], o);
        for (let e = 0; e < a.members.length; e++)
          (n.add(a.members[e]), r.set(a.members[e], 0));
      }
      let i = 0,
        a = 0;
      for (let o = 0; o < this.graph.edges.length; o++) {
        let s = this.graph.edges[o];
        if (!n.has(s.sourceId) || !n.has(s.targetId)) continue;
        let c = 0;
        for (let n = 0; n < e.length; n++) {
          let r = t.get(e[n]);
          r.has(s.sourceId) && r.has(s.targetId) && c++;
        }
        c !== 0 &&
          (i++,
          r.set(s.sourceId, r.get(s.sourceId) + 1),
          r.set(s.targetId, r.get(s.targetId) + 1),
          c === 1 && a++);
      }
      let o = 0;
      for (let e of r.values()) e !== 3 && o++;
      return {
        atomCount: n.size,
        edgeCount: i,
        nonCageAtomCount: o,
        boundaryEdgeRatio: i === 0 ? 1 : a / i,
      };
    }
    createBridgedRing(e, t) {
      let n = new Set(),
        r = new Set(),
        i = new Set();
      for (let t = 0; t < e.length; t++) {
        let n = this.getRing(e[t]);
        n.isPartOfBridged = !0;
        for (let e = 0; e < n.members.length; e++) r.add(n.members[e]);
        for (let t = 0; t < n.neighbours.length; t++) {
          let r = n.neighbours[t];
          e.indexOf(r) === -1 && i.add(n.neighbours[t]);
        }
      }
      let a = new Set();
      for (let t of r) {
        let r = this.graph.vertices[t],
          i = f.intersection(e, r.value.rings);
        r.value.rings.length === 1 || i.length === 1
          ? n.add(r.id)
          : a.add(r.id);
      }
      let o = [];
      for (let e of a) {
        let t = this.graph.vertices[e],
          r = !1;
        for (let e = 0; e < t.edges.length; e++)
          this.edgeRingCount(t.edges[e]) === 1 && (r = !0);
        r
          ? ((t.value.isBridgeNode = !0), n.add(t.id))
          : ((t.value.isBridge = !0), o.push(t.id), n.add(t.id));
      }
      let s = new b([...n]);
      (this.addRing(s),
        (s.isBridged = !0),
        (s.insiders = o),
        (s.neighbours = [...i]));
      for (let t = 0; t < e.length; t++)
        s.rings.push(this.getRing(e[t]).clone());
      for (let e = 0; e < s.members.length; e++)
        this.graph.vertices[s.members[e]].value.bridgedRing = s.id;
      for (let e = 0; e < o.length; e++) {
        let t = this.graph.vertices[o[e]];
        t.value.rings = [];
      }
      for (let t of n) {
        let n = this.graph.vertices[t];
        ((n.value.rings = f.removeAll(n.value.rings, e)),
          n.value.rings.push(s.id));
      }
      for (let t = 0; t < e.length; t++)
        for (let n = t + 1; n < e.length; n++)
          this.removeRingConnectionsBetween(e[t], e[n]);
      for (let t of i) {
        let n = this.getRingConnections(t, e);
        for (let e = 0; e < n.length; e++)
          this.getRingConnection(n[e]).updateOther(s.id, t);
        this.getRing(t).neighbours.push(s.id);
      }
      return s;
    }
    areVerticesInSameRing(e, t) {
      for (let n = 0; n < e.value.rings.length; n++)
        for (let r = 0; r < t.value.rings.length; r++)
          if (e.value.rings[n] === t.value.rings[r]) return !0;
      return !1;
    }
    getCommonRings(e, t) {
      let n = [];
      for (let r = 0; r < e.value.rings.length; r++)
        for (let i = 0; i < t.value.rings.length; i++)
          e.value.rings[r] == t.value.rings[i] && n.push(e.value.rings[r]);
      return n;
    }
    getLargestOrAromaticCommonRing(e, t) {
      let n = this.getCommonRings(e, t),
        r = 0,
        i = null;
      for (let e = 0; e < n.length; e++) {
        let t = this.getRing(n[e]),
          a = t.getSize();
        if (t.isBenzeneLike(this.graph.vertices)) return t;
        a > r && ((r = a), (i = t));
      }
      return i;
    }
    getVerticesAt(e, t, n) {
      let r = [];
      for (let i = 0; i < this.graph.vertices.length; i++) {
        let a = this.graph.vertices[i];
        a.id === n ||
          !a.positioned ||
          (e.distanceSq(a.position) <= t * t && r.push(a.id));
      }
      return r;
    }
    getClosestVertex(e) {
      let t = 99999,
        n = null;
      for (let r = 0; r < this.graph.vertices.length; r++) {
        let i = this.graph.vertices[r];
        if (i.id === e.id) continue;
        let a = e.position.distanceSq(i.position);
        a < t && ((t = a), (n = i));
      }
      return n;
    }
    addRing(e) {
      return ((e.id = this.ringIdCounter++), this.rings.push(e), e.id);
    }
    removeRing(e) {
      ((this.rings = this.rings.filter(function (t) {
        return t.id !== e;
      })),
        (this.ringConnections = this.ringConnections.filter(function (t) {
          return t.firstRingId !== e && t.secondRingId !== e;
        })));
      for (let t = 0; t < this.rings.length; t++) {
        let n = this.rings[t];
        n.neighbours = n.neighbours.filter(function (t) {
          return t !== e;
        });
      }
    }
    getRing(e) {
      for (let t = 0; t < this.rings.length; t++)
        if (this.rings[t].id == e) return this.rings[t];
    }
    addRingConnection(e) {
      return (
        (e.id = this.ringConnectionIdCounter++),
        this.ringConnections.push(e),
        e.id
      );
    }
    removeRingConnection(e) {
      this.ringConnections = this.ringConnections.filter(function (t) {
        return t.id !== e;
      });
    }
    removeRingConnectionsBetween(e, t) {
      let n = [];
      for (let r = 0; r < this.ringConnections.length; r++) {
        let i = this.ringConnections[r];
        ((i.firstRingId === e && i.secondRingId === t) ||
          (i.firstRingId === t && i.secondRingId === e)) &&
          n.push(i.id);
      }
      for (let e = 0; e < n.length; e++) this.removeRingConnection(n[e]);
    }
    getRingConnection(e) {
      for (let t = 0; t < this.ringConnections.length; t++)
        if (this.ringConnections[t].id == e) return this.ringConnections[t];
    }
    getRingConnections(e, t) {
      let n = [];
      for (let r = 0; r < this.ringConnections.length; r++) {
        let i = this.ringConnections[r];
        for (let r = 0; r < t.length; r++) {
          let a = t[r];
          ((i.firstRingId === e && i.secondRingId === a) ||
            (i.firstRingId === a && i.secondRingId === e)) &&
            n.push(i.id);
        }
      }
      return n;
    }
    getOverlapScore() {
      let e = 0,
        t = new Float32Array(this.graph.vertices.length);
      for (let e = 0; e < this.graph.vertices.length; e++) t[e] = 0;
      for (let n = 0; n < this.graph.vertices.length; n++) {
        let r = this.graph.vertices.length;
        for (; --r > n;) {
          let i = this.graph.vertices[n],
            a = this.graph.vertices[r];
          if (!i.value.isDrawn || !a.value.isDrawn) continue;
          let o = h.subtract(i.position, a.position).lengthSq();
          if (o < this.opts.bondLengthSq) {
            let i =
              (this.opts.bondLength - Math.sqrt(o)) / this.opts.bondLength;
            ((e += i), (t[n] += i), (t[r] += i));
          }
        }
      }
      let n = [];
      for (let e = 0; e < this.graph.vertices.length; e++)
        n.push({ id: e, score: t[e] });
      return (
        n.sort(function (e, t) {
          return t.score - e.score;
        }),
        { total: e, scores: n, vertexScores: t }
      );
    }
    chooseSide(e, t, n) {
      let r = e.getNeighbours(t.id),
        i = t.getNeighbours(e.id),
        a = r.length,
        o = i.length,
        s = f.merge(r, i),
        c = [0, 0];
      for (let r = 0; r < s.length; r++)
        this.graph.vertices[s[r]].position.sameSideAs(
          e.position,
          t.position,
          n[0],
        )
          ? c[0]++
          : c[1]++;
      let l = [0, 0];
      for (let r = 0; r < this.graph.vertices.length; r++)
        this.graph.vertices[r].position.sameSideAs(e.position, t.position, n[0])
          ? l[0]++
          : l[1]++;
      return {
        totalSideCount: l,
        totalPosition: l[0] > l[1] ? 0 : 1,
        sideCount: c,
        position: c[0] > c[1] ? 0 : 1,
        anCount: a,
        bnCount: o,
      };
    }
    setRingCenter(e) {
      let t = e.getSize(),
        n = new h(0, 0);
      for (let r = 0; r < t; r++)
        n.add(this.graph.vertices[e.members[r]].position);
      e.center = n.divide(t);
    }
    getSubringCenter(e, t) {
      let n = t.value.originalRings,
        r = e.center,
        i = Number.MAX_VALUE;
      for (let t = 0; t < n.length; t++)
        for (let a = 0; a < e.rings.length; a++)
          n[t] === e.rings[a].id &&
            e.rings[a].getSize() < i &&
            ((r = e.rings[a].center), (i = e.rings[a].getSize()));
      return r;
    }
    drawEdges(e) {
      let t = Array(this.graph.edges.length);
      (t.fill(!1),
        this.graph.traverseBF(0, (n) => {
          let r = this.graph.getEdges(n.id);
          for (let n = 0; n < r.length; n++) {
            let i = r[n];
            t[i] || ((t[i] = !0), this.drawEdge(i, e));
          }
        }));
      for (let e = 0; e < this.rings.length; e++) {
        let t = this.rings[e];
        this.isRingAromatic(t) &&
          ((t.isPartOfBridged && !this.isRingRegularPolygon(t)) ||
            this.canvasWrapper.drawAromaticityRing(t));
      }
    }
    drawEdge(e, t) {
      let n = this.graph.edges[e],
        r = this.graph.vertices[n.sourceId],
        i = this.graph.vertices[n.targetId],
        a = r.value.element,
        o = i.value.element;
      if (
        (!r.value.isDrawn || !i.value.isDrawn) &&
        this.opts.atomVisualization === "default"
      )
        return;
      let s = r.position,
        c = i.position,
        l = this.getEdgeNormals(n),
        u = f.clone(l);
      (u[0].multiplyScalar(10).add(s), u[1].multiplyScalar(10).add(s));
      let d =
        n.isPartOfAromaticRing &&
        r.value.bridgedRing !== null &&
        i.value.bridgedRing !== null
          ? this.getLargestOrAromaticCommonRing(r, i)
          : null;
      if (
        n.bondType === `=` ||
        this.getRingbondType(r, i) === `=` ||
        (d && !this.isRingRegularPolygon(d))
      ) {
        let e = this.areVerticesInSameRing(r, i),
          t = this.chooseSide(r, i, u);
        if (e) {
          let e = this.getLargestOrAromaticCommonRing(r, i).center;
          (l[0].multiplyScalar(this.opts.bondSpacing),
            l[1].multiplyScalar(this.opts.bondSpacing));
          let t = null;
          ((t = e.sameSideAs(r.position, i.position, h.add(s, l[0]))
            ? new g(h.add(s, l[0]), h.add(c, l[0]), a, o)
            : new g(h.add(s, l[1]), h.add(c, l[1]), a, o)),
            t.shorten(
              this.opts.bondLength -
                this.opts.shortBondLength * this.opts.bondLength,
            ),
            n.isPartOfAromaticRing
              ? this.canvasWrapper.drawLine(t, !0)
              : this.canvasWrapper.drawLine(t),
            this.canvasWrapper.drawLine(new g(s, c, a, o)));
        } else if (n.center || (r.isTerminal() && i.isTerminal())) {
          (l[0].multiplyScalar(this.opts.halfBondSpacing),
            l[1].multiplyScalar(this.opts.halfBondSpacing));
          let e = new g(h.add(s, l[0]), h.add(c, l[0]), a, o),
            t = new g(h.add(s, l[1]), h.add(c, l[1]), a, o);
          (this.canvasWrapper.drawLine(e), this.canvasWrapper.drawLine(t));
        } else if (
          (t.anCount == 0 && t.bnCount > 1) ||
          (t.bnCount == 0 && t.anCount > 1)
        ) {
          (l[0].multiplyScalar(this.opts.halfBondSpacing),
            l[1].multiplyScalar(this.opts.halfBondSpacing));
          let e = new g(h.add(s, l[0]), h.add(c, l[0]), a, o),
            t = new g(h.add(s, l[1]), h.add(c, l[1]), a, o);
          (this.canvasWrapper.drawLine(e), this.canvasWrapper.drawLine(t));
        } else if (t.sideCount[0] > t.sideCount[1]) {
          (l[0].multiplyScalar(this.opts.bondSpacing),
            l[1].multiplyScalar(this.opts.bondSpacing));
          let e = new g(h.add(s, l[0]), h.add(c, l[0]), a, o);
          (e.shorten(
            this.opts.bondLength -
              this.opts.shortBondLength * this.opts.bondLength,
          ),
            this.canvasWrapper.drawLine(e),
            this.canvasWrapper.drawLine(new g(s, c, a, o)));
        } else if (t.sideCount[0] < t.sideCount[1]) {
          (l[0].multiplyScalar(this.opts.bondSpacing),
            l[1].multiplyScalar(this.opts.bondSpacing));
          let e = new g(h.add(s, l[1]), h.add(c, l[1]), a, o);
          (e.shorten(
            this.opts.bondLength -
              this.opts.shortBondLength * this.opts.bondLength,
          ),
            this.canvasWrapper.drawLine(e),
            this.canvasWrapper.drawLine(new g(s, c, a, o)));
        } else if (t.totalSideCount[0] > t.totalSideCount[1]) {
          (l[0].multiplyScalar(this.opts.bondSpacing),
            l[1].multiplyScalar(this.opts.bondSpacing));
          let e = new g(h.add(s, l[0]), h.add(c, l[0]), a, o);
          (e.shorten(
            this.opts.bondLength -
              this.opts.shortBondLength * this.opts.bondLength,
          ),
            this.canvasWrapper.drawLine(e),
            this.canvasWrapper.drawLine(new g(s, c, a, o)));
        } else if (t.totalSideCount[0] <= t.totalSideCount[1]) {
          (l[0].multiplyScalar(this.opts.bondSpacing),
            l[1].multiplyScalar(this.opts.bondSpacing));
          let e = new g(h.add(s, l[1]), h.add(c, l[1]), a, o);
          (e.shorten(
            this.opts.bondLength -
              this.opts.shortBondLength * this.opts.bondLength,
          ),
            this.canvasWrapper.drawLine(e),
            this.canvasWrapper.drawLine(new g(s, c, a, o)));
        }
      } else if (n.bondType === `#`) {
        (l[0].multiplyScalar(this.opts.bondSpacing / 1.5),
          l[1].multiplyScalar(this.opts.bondSpacing / 1.5));
        let e = new g(h.add(s, l[0]), h.add(c, l[0]), a, o),
          t = new g(h.add(s, l[1]), h.add(c, l[1]), a, o);
        (this.canvasWrapper.drawLine(e),
          this.canvasWrapper.drawLine(t),
          this.canvasWrapper.drawLine(new g(s, c, a, o)));
      } else if (n.bondType !== `.`) {
        let e = r.value.isStereoCenter,
          t = i.value.isStereoCenter;
        n.wedge === `up`
          ? this.canvasWrapper.drawWedge(new g(s, c, a, o, e, t))
          : n.wedge === `down`
            ? this.canvasWrapper.drawDashedWedge(new g(s, c, a, o, e, t))
            : this.canvasWrapper.drawLine(new g(s, c, a, o, e, t));
      }
      if (t) {
        let t = h.midpoint(s, c);
        this.canvasWrapper.drawDebugText(t.x, t.y, `e: ` + e);
      }
    }
    drawVertices(t) {
      for (let n = 0; n < this.graph.vertices.length; n++) {
        let r = this.graph.vertices[n],
          i = r.value,
          a = 0,
          o = 0,
          s = i.element,
          c = i.countImplicitHydrogens(),
          l = r.getTextDirection(this.graph.vertices),
          u = e.getEffectiveShowCarbonsMode(this.opts),
          d =
            u === `terminal` || s !== `C` || i.hasAttachedPseudoElements
              ? r.isTerminal()
              : !1,
          p = i.element === `C`;
        if (s === `C`) {
          let e = i.rings && i.rings.length > 0;
          u === `none`
            ? ((p = !0), (d = !1))
            : (u === `all` || (u === `acyclic` && !e)) && ((p = !1), (d = !0));
        }
        if (
          (i.element === `N` && i.isPartOfAromaticRing && (c = 0),
          i.bracket &&
            ((c = i.bracket.hcount),
            (a = i.bracket.charge),
            (o = i.bracket.isotope)),
          (a || o || this.graph.vertices.length < 3) && (p = !1),
          this.opts.atomVisualization === `allballs`)
        )
          this.canvasWrapper.drawBall(r.position.x, r.position.y, s);
        else if (
          (i.isDrawn &&
            (!p || i.drawExplicit || d || i.hasAttachedPseudoElements)) ||
          this.graph.vertices.length === 1
        )
          this.opts.atomVisualization === "default"
            ? this.canvasWrapper.drawText(
                r.position.x,
                r.position.y,
                s,
                c,
                l,
                d,
                a,
                o,
                this.graph.vertices.length,
                i.getAttachedPseudoElements(),
              )
            : this.opts.atomVisualization === `balls` &&
              this.canvasWrapper.drawBall(r.position.x, r.position.y, s);
        else if (r.getNeighbourCount() === 2 && r.forcePositioned == 1) {
          let e = this.graph.vertices[r.neighbours[0]].position,
            t = this.graph.vertices[r.neighbours[1]].position,
            n = h.threePointangle(r.position, e, t);
          Math.abs(Math.PI - n) < 0.1 &&
            this.canvasWrapper.drawPoint(r.position.x, r.position.y, s);
        }
        if (t) {
          let e = `v: ` + r.id + ` ` + f.print(i.ringbonds);
          this.canvasWrapper.drawDebugText(r.position.x, r.position.y, e);
        }
      }
      if (this.opts.debug)
        for (let e = 0; e < this.rings.length; e++) {
          let t = this.rings[e].center;
          this.canvasWrapper.drawDebugPoint(t.x, t.y, `r: ` + this.rings[e].id);
        }
    }
    position() {
      let e = null;
      for (let t = 0; t < this.graph.vertices.length; t++)
        if (this.graph.vertices[t].value.bridgedRing !== null) {
          e = this.graph.vertices[t];
          break;
        }
      for (let t = 0; t < this.rings.length; t++)
        this.rings[t].isBridged &&
          (e = this.graph.vertices[this.rings[t].members[0]]);
      (this.rings.length > 0 &&
        e === null &&
        (e = this.graph.vertices[this.rings[0].members[0]]),
        e === null && (e = this.graph.vertices[0]),
        this.createNextBond(e, null, 0));
    }
    backupRingInformation() {
      ((this.originalRings = []), (this.originalRingConnections = []));
      for (let e = 0; e < this.rings.length; e++)
        this.originalRings.push(this.rings[e]);
      for (let e = 0; e < this.ringConnections.length; e++)
        this.originalRingConnections.push(this.ringConnections[e]);
      for (let e = 0; e < this.graph.vertices.length; e++)
        this.graph.vertices[e].value.backupRings();
    }
    restoreRingInformation() {
      let e = this.getBridgedRings();
      ((this.rings = []), (this.ringConnections = []));
      for (let t = 0; t < e.length; t++) {
        let n = e[t];
        for (let e = 0; e < n.rings.length; e++) {
          let t = n.rings[e];
          this.originalRings[t.id].center = t.center;
        }
      }
      for (let e = 0; e < this.originalRings.length; e++)
        this.rings.push(this.originalRings[e]);
      for (let e = 0; e < this.originalRingConnections.length; e++)
        this.ringConnections.push(this.originalRingConnections[e]);
      for (let e = 0; e < this.graph.vertices.length; e++)
        this.graph.vertices[e].value.restoreRings();
    }
    createRing(e, t = null, n = null, r = null) {
      if (e.positioned) return;
      t ||= new h(0, 0);
      let i = e.getOrderedNeighbours(this.ringConnections),
        a = n ? h.subtract(n.position, t).angle() : 0,
        o = _.polyCircumradius(this.opts.bondLength, e.getSize()),
        s = _.centralAngle(e.getSize());
      e.centralAngle = s;
      let c = a,
        l = n ? n.id : null;
      if (
        (e.members.indexOf(l) === -1 &&
          (n && (n.positioned = !1), (l = e.members[0])),
        e.isBridged)
      ) {
        (this.graph.kkLayout(
          e.members.slice(),
          t,
          n.id,
          e,
          this.opts.bondLength,
          this.opts.kkThreshold,
          this.opts.kkInnerThreshold,
          this.opts.kkMaxIteration,
          this.opts.kkMaxInnerIteration,
          this.opts.kkMaxEnergy,
        ),
          (e.positioned = !0),
          this.setRingCenter(e),
          (t = e.center));
        for (let t = 0; t < e.rings.length; t++) this.setRingCenter(e.rings[t]);
      } else
        e.eachMember(
          this.graph.vertices,
          (n) => {
            let r = this.graph.vertices[n];
            (r.positioned ||
              r.setPosition(t.x + Math.cos(c) * o, t.y + Math.sin(c) * o),
              (c += s),
              (!e.isBridged || e.rings.length < 3) &&
                ((r.angle = c), (r.positioned = !0)));
          },
          l,
          r ? r.id : null,
        );
      ((e.positioned = !0), (e.center = t));
      for (let n = 0; n < i.length; n++) {
        let r = this.getRing(i[n].neighbour);
        if (r.positioned) continue;
        let a = y.getVertices(this.ringConnections, e.id, r.id);
        if (a.length === 2) {
          ((e.isFused = !0), (r.isFused = !0));
          let n = this.graph.vertices[a[0]],
            i = this.graph.vertices[a[1]],
            o = h.midpoint(n.position, i.position),
            s = h.normals(n.position, i.position);
          (s[0].normalize(), s[1].normalize());
          let c = _.polyCircumradius(this.opts.bondLength, r.getSize()),
            l = _.apothem(c, r.getSize());
          (s[0].multiplyScalar(l).add(o), s[1].multiplyScalar(l).add(o));
          let u = s[0];
          h.subtract(t, s[1]).lengthSq() > h.subtract(t, s[0]).lengthSq() &&
            (u = s[1]);
          let d = h.subtract(n.position, u),
            f = h.subtract(i.position, u);
          d.clockwise(f) === -1
            ? r.positioned || this.createRing(r, u, n, i)
            : r.positioned || this.createRing(r, u, i, n);
        } else if (a.length === 1) {
          ((e.isSpiro = !0), (r.isSpiro = !0));
          let n = this.graph.vertices[a[0]],
            i = h.subtract(t, n.position);
          (i.invert(), i.normalize());
          let o = _.polyCircumradius(this.opts.bondLength, r.getSize());
          (i.multiplyScalar(o),
            i.add(n.position),
            r.positioned || this.createRing(r, i, n));
        }
      }
      for (let t = 0; t < e.members.length; t++) {
        let n = this.graph.vertices[e.members[t]],
          r = n.neighbours;
        for (let e = 0; e < r.length; e++) {
          let t = this.graph.vertices[r[e]];
          t.positioned ||
            ((t.value.isConnectedToRing = !0), this.createNextBond(t, n, 0));
        }
      }
    }
    fixDoubleBondStereo() {
      let e = this.graph;
      for (let t = 0; t < e.edges.length; t++) {
        let n = e.edges[t];
        if (n.bondType !== `=`) continue;
        let r = n.sourceId,
          i = n.targetId;
        if (this.areVerticesInSameRing(e.vertices[r], e.vertices[i])) continue;
        let a = null,
          o = null;
        for (let t of e.vertices[r].getNeighbours()) {
          if (t === i) continue;
          let n = e.getEdge(r, t);
          if (n && (n.bondType === `/` || n.bondType === `\\`)) {
            a = {
              nid: t,
              above:
                n.sourceId === r ? n.bondType === `/` : n.bondType === `\\`,
            };
            break;
          }
        }
        for (let t of e.vertices[i].getNeighbours()) {
          if (t === r) continue;
          let n = e.getEdge(i, t);
          if (n && (n.bondType === `/` || n.bondType === `\\`)) {
            o = {
              nid: t,
              above:
                n.sourceId === i ? n.bondType === `/` : n.bondType === `\\`,
            };
            break;
          }
        }
        if (!a || !o) continue;
        let s = a.above === o.above,
          c = e.vertices[r].position,
          l = e.vertices[i].position,
          u = e.vertices[a.nid].position,
          d = e.vertices[o.nid].position,
          f = l.x - c.x,
          p = l.y - c.y,
          m = f * (u.y - c.y) - p * (u.x - c.x),
          h = f * (d.y - l.y) - p * (d.x - l.x);
        if (s === (m > 0 == h > 0)) continue;
        let g = 0,
          _ = 0;
        for (let t of e.vertices[r].getNeighbours()) {
          if (t === i) continue;
          let n = e.getEdge(r, t);
          n && (n.bondType === `/` || n.bondType === `\\`) && g++;
        }
        for (let t of e.vertices[i].getNeighbours()) {
          if (t === r) continue;
          let n = e.getEdge(i, t);
          n && (n.bondType === `/` || n.bondType === `\\`) && _++;
        }
        let v = g <= _ ? r : i,
          y = g <= _ ? i : r,
          b = e.vertices[y].position,
          x = f * f + p * p;
        x < 0.001 ||
          e.traverseTree(v, y, (e) => {
            let t = e.position.x - b.x,
              n = e.position.y - b.y,
              r = t * f + n * p;
            ((e.position.x = b.x + (2 * r * f) / x - t),
              (e.position.y = b.y + (2 * r * p) / x - n));
            for (let t = 0; t < e.value.anchoredRings.length; t++) {
              let n = this.rings[e.value.anchoredRings[t]];
              if (n) {
                let e = n.center.x - b.x,
                  t = n.center.y - b.y,
                  r = e * f + t * p;
                ((n.center.x = b.x + (2 * r * f) / x - e),
                  (n.center.y = b.y + (2 * r * p) / x - t));
              }
            }
          });
      }
    }
    rotateSubtree(e, t, n, r) {
      let i = r.clone();
      this.graph.traverseTree(e, t, (e) => {
        e.position.rotateAround(n, i);
        for (let t = 0; t < e.value.anchoredRings.length; t++) {
          let r = this.rings[e.value.anchoredRings[t]];
          r && r.center.rotateAround(n, i);
        }
      });
    }
    getSubtreeOverlapScore(e, t, n) {
      let r = 0,
        i = new h(0, 0),
        a = 0;
      return (
        this.graph.traverseTree(e, t, (e) => {
          if (!e.value.isDrawn) return;
          let t = n[e.id];
          t > this.opts.overlapSensitivity && ((r += t), a++);
          let o = this.graph.vertices[e.id].position.clone();
          (o.multiplyScalar(t), i.add(o));
        }),
        i.divide(r),
        { value: r / a, center: i }
      );
    }
    getCurrentCenterOfMass() {
      let e = new h(0, 0),
        t = 0;
      for (let n = 0; n < this.graph.vertices.length; n++) {
        let r = this.graph.vertices[n];
        r.positioned && (e.add(r.position), t++);
      }
      return e.divide(t);
    }
    getCurrentCenterOfMassInNeigbourhood(e, t = this.opts.bondLength * 2) {
      let n = new h(0, 0),
        r = 0,
        i = t * t;
      for (let t = 0; t < this.graph.vertices.length; t++) {
        let a = this.graph.vertices[t];
        a.positioned &&
          e.distanceSq(a.position) < i &&
          (n.add(a.position), r++);
      }
      return n.divide(r);
    }
    resolvePrimaryOverlaps() {
      let e = [],
        t = Array(this.graph.vertices.length);
      for (let n = 0; n < this.rings.length; n++) {
        let r = this.rings[n];
        for (let n = 0; n < r.members.length; n++) {
          let i = this.graph.vertices[r.members[n]];
          if (t[i.id]) continue;
          t[i.id] = !0;
          let a = this.getNonRingNeighbours(i.id);
          if (a.length > 1) {
            let t = [];
            for (let e = 0; e < i.value.rings.length; e++)
              t.push(i.value.rings[e]);
            e.push({ common: i, rings: t, vertices: a });
          } else if (a.length === 1 && i.value.rings.length === 2) {
            let t = [];
            for (let e = 0; e < i.value.rings.length; e++)
              t.push(i.value.rings[e]);
            e.push({ common: i, rings: t, vertices: a });
          }
        }
      }
      for (let t = 0; t < e.length; t++) {
        let n = e[t];
        if (n.vertices.length === 2) {
          let e = n.vertices[0],
            t = n.vertices[1];
          if (!e.value.isDrawn || !t.value.isDrawn) continue;
          let r = (2 * Math.PI - this.getRing(n.rings[0]).getAngle()) / 6;
          (this.rotateSubtree(e.id, n.common.id, r, n.common.position),
            this.rotateSubtree(t.id, n.common.id, -r, n.common.position));
          let i = this.getOverlapScore(),
            a = this.getSubtreeOverlapScore(e.id, n.common.id, i.vertexScores),
            o = this.getSubtreeOverlapScore(t.id, n.common.id, i.vertexScores),
            s = a.value + o.value;
          (this.rotateSubtree(e.id, n.common.id, -2 * r, n.common.position),
            this.rotateSubtree(t.id, n.common.id, 2 * r, n.common.position),
            (i = this.getOverlapScore()),
            (a = this.getSubtreeOverlapScore(
              e.id,
              n.common.id,
              i.vertexScores,
            )),
            (o = this.getSubtreeOverlapScore(
              t.id,
              n.common.id,
              i.vertexScores,
            )),
            a.value + o.value > s &&
              (this.rotateSubtree(e.id, n.common.id, 2 * r, n.common.position),
              this.rotateSubtree(
                t.id,
                n.common.id,
                -2 * r,
                n.common.position,
              )));
        } else n.vertices.length === 1 && n.rings.length;
      }
    }
    resolveSecondaryOverlaps(e) {
      for (let t = 0; t < e.length; t++)
        if (e[t].score > this.opts.overlapSensitivity) {
          let n = this.graph.vertices[e[t].id];
          if (n.isTerminal()) {
            let e = this.getClosestVertex(n);
            if (e) {
              let t = null;
              t = e.isTerminal()
                ? e.id === 0
                  ? this.graph.vertices[1].position
                  : e.previousPosition
                : e.id === 0
                  ? this.graph.vertices[1].position
                  : e.position;
              let r =
                n.id === 0
                  ? this.graph.vertices[1].position
                  : n.previousPosition;
              n.position.rotateAwayFrom(t, r, _.toRad(20));
            }
          }
        }
    }
    getLastAngle(e) {
      for (; e;) {
        let t = this.graph.vertices[e];
        if (t.value.rings.length > 0) return 0;
        if (t.angle) return t.angle;
        e = t.parentVertexId;
      }
      return 0;
    }
    createNextBond(e, t = null, n = 0, r = !1, i = !1) {
      if (e.positioned && !i) return;
      let a = !1;
      if (t) {
        let n = this.graph.getEdge(e.id, t.id);
        (n.bondType === `/` || n.bondType === `\\`) &&
          ++this.doubleBondConfigCount % 2 == 1 &&
          this.doubleBondConfig === null &&
          ((this.doubleBondConfig = n.bondType),
          (a = !0),
          t.parentVertexId === null &&
            e.value.branchBond &&
            (this.doubleBondConfig === `/`
              ? (this.doubleBondConfig = `\\`)
              : this.doubleBondConfig === `\\` &&
                (this.doubleBondConfig = `/`)));
      }
      if (!i) {
        if (t) {
          if (t.value.rings.length > 0) {
            let n = t.neighbours,
              r = null,
              i = new h(0, 0);
            if (t.value.bridgedRing === null && t.value.rings.length > 1)
              for (let e = 0; e < n.length; e++) {
                let i = this.graph.vertices[n[e]];
                if (f.containsAll(i.value.rings, t.value.rings)) {
                  r = i;
                  break;
                }
              }
            if (r === null) {
              for (let e = 0; e < n.length; e++) {
                let r = this.graph.vertices[n[e]];
                r.positioned &&
                  this.areVerticesInSameRing(r, t) &&
                  i.add(h.subtract(r.position, t.position));
              }
              if (i.lengthSq() < 1) {
                let e = null;
                ((e =
                  t.value.bridgedRing === null
                    ? this.getRing(t.value.rings[0])
                    : this.getRing(t.value.bridgedRing)),
                  (i =
                    e && e.center
                      ? h.subtract(e.center, t.position)
                      : new h(1, 0)));
              }
              i.invert()
                .normalize()
                .multiplyScalar(this.opts.bondLength)
                .add(t.position);
            } else i = r.position.clone().rotateAround(Math.PI, t.position);
            ((e.previousPosition = t.position),
              e.setPositionFromVector(i),
              (e.positioned = !0));
          } else {
            let r = new h(this.opts.bondLength, 0);
            (r.rotate(n),
              r.add(t.position),
              e.setPositionFromVector(r),
              (e.previousPosition = t.position),
              (e.positioned = !0));
          }
        } else {
          let t = new h(this.opts.bondLength, 0);
          (t.rotate(_.toRad(-60)),
            (e.previousPosition = t),
            e.setPosition(this.opts.bondLength, 0),
            (e.angle = _.toRad(-60)),
            e.value.bridgedRing === null && (e.positioned = !0));
        }
      }
      if (e.value.bridgedRing !== null) {
        let t = this.getRing(e.value.bridgedRing);
        if (!t.positioned) {
          let n = h.subtract(e.previousPosition, e.position);
          (n.invert(), n.normalize());
          let r = _.polyCircumradius(this.opts.bondLength, t.members.length);
          (n.multiplyScalar(r), n.add(e.position), this.createRing(t, n, e));
        }
      } else if (e.value.rings.length > 0) {
        let t = this.getRing(e.value.rings[0]);
        if (!t.positioned) {
          let n = h.subtract(e.previousPosition, e.position);
          (n.invert(), n.normalize());
          let r = _.polyCircumradius(this.opts.bondLength, t.getSize());
          (n.multiplyScalar(r), n.add(e.position), this.createRing(t, n, e));
        }
      } else {
        let n = e.getNeighbours(),
          i = [];
        for (let e = 0; e < n.length; e++)
          this.graph.vertices[n[e]].value.isDrawn && i.push(n[e]);
        t && (i = f.remove(i, t.id));
        let o = e.getAngle();
        if (i.length === 1) {
          let n = this.graph.vertices[i[0]],
            s = t ? this.graph.getEdge(e.id, t.id) : null,
            c = this.graph.getEdge(e.id, n.id);
          if (s && c && s.weight + c.weight >= 4)
            ((s.center = !0),
              (c.center = !0),
              (n.angle = 0),
              s.weight === c.weight && (e.value.drawExplicit = !0),
              this.createNextBond(n, e, o + n.angle));
          else if (t && t.value.rings.length > 0) {
            let t = _.toRad(60),
              r = -t,
              i = new h(this.opts.bondLength, 0),
              a = new h(this.opts.bondLength, 0);
            (i.rotate(t).add(e.position), a.rotate(r).add(e.position));
            let s = this.getCurrentCenterOfMass();
            ((n.angle = i.distanceSq(s) < a.distanceSq(s) ? r : t),
              this.createNextBond(n, e, o + n.angle));
          } else {
            let i = this.getLastAngle(e.id);
            if (((i = i >= 0 ? 1.0472 : -1.0472), t && !a)) {
              let t = this.graph.getEdge(e.id, n.id).bondType;
              t === `/`
                ? (this.doubleBondConfig === `/` ||
                    (this.doubleBondConfig === `\\` && (i = -i)),
                  (this.doubleBondConfig = null))
                : t === `\\` &&
                  (this.doubleBondConfig === `/`
                    ? (i = -i)
                    : this.doubleBondConfig,
                  (this.doubleBondConfig = null));
            }
            ((n.angle = r ? i : -i), this.createNextBond(n, e, o + n.angle));
          }
        } else if (i.length === 2) {
          let n = e.angle;
          n ||= 1.0472;
          let r = this.graph.getTreeDepth(i[0], e.id),
            a = this.graph.getTreeDepth(i[1], e.id),
            s = this.graph.vertices[i[0]],
            c = this.graph.vertices[i[1]];
          ((s.value.subtreeDepth = r), (c.value.subtreeDepth = a));
          let l = this.graph.getTreeDepth(t ? t.id : null, e.id);
          t && (t.value.subtreeDepth = l);
          let u = 0,
            d = 1;
          c.value.element === `C` && s.value.element !== `C` && a > 1 && r < 5
            ? ((u = 1), (d = 0))
            : c.value.element !== `C` &&
                s.value.element === `C` &&
                r > 1 &&
                a < 5
              ? ((u = 0), (d = 1))
              : a > r && ((u = 1), (d = 0));
          let f = this.graph.vertices[i[u]],
            p = this.graph.vertices[i[d]],
            m = l < r && l < a;
          ((p.angle = n),
            (f.angle = -n),
            this.doubleBondConfig === `\\`
              ? p.value.branchBond === `\\` && ((p.angle = -n), (f.angle = n))
              : this.doubleBondConfig === `/` &&
                p.value.branchBond === `/` &&
                ((p.angle = -n), (f.angle = n)),
            this.createNextBond(p, e, o + p.angle, m),
            this.createNextBond(f, e, o + f.angle, m));
        } else if (i.length > 0) {
          let n = i.map((t) => {
            let n = this.graph.vertices[t],
              r = this.graph.getTreeDepth(t, e.id);
            return ((n.value.subtreeDepth = r), n);
          });
          if (
            (n.sort((e, t) => t.value.subtreeDepth - e.value.subtreeDepth),
            i.length === 3 &&
              t &&
              t.parentVertexId !== null &&
              t.value.rings.length < 1 &&
              n[2].value.rings.length < 1 &&
              n[1].value.rings.length < 1 &&
              n[0].value.rings.length < 1 &&
              n[2].value.subtreeDepth === 1 &&
              n[1].value.subtreeDepth === 1 &&
              n[0].value.subtreeDepth > 1)
          )
            (e.angle >= 0
              ? ((n[0].angle = -1.0472),
                (n[1].angle = _.toRad(30)),
                (n[2].angle = _.toRad(90)))
              : ((n[0].angle = 1.0472),
                (n[1].angle = -_.toRad(30)),
                (n[2].angle = -_.toRad(90))),
              this.createNextBond(n[0], e, o + n[0].angle),
              this.createNextBond(n[1], e, o + n[1].angle),
              this.createNextBond(n[2], e, o + n[2].angle));
          else {
            let r = i.length + +!!t,
              a = (2 * Math.PI) / r,
              s = a,
              c = 0;
            for (
              i.length % 2 == 0
                ? (s /= 2)
                : (this.createNextBond(n[0], e, o), (c = 1));
              c < i.length;
            )
              (this.createNextBond(n[c + 0], e, o + s),
                this.createNextBond(n[c + 1], e, o - s),
                (s += a),
                (c += 2));
          }
        }
      }
    }
    getCommonRingbondNeighbour(e) {
      let t = e.neighbours;
      for (let n = 0; n < t.length; n++) {
        let r = this.graph.vertices[t[n]];
        if (f.containsAll(r.value.rings, e.value.rings)) return r;
      }
      return null;
    }
    isPointInRing(e) {
      for (let t = 0; t < this.rings.length; t++) {
        let n = this.rings[t];
        if (!n.positioned) continue;
        let r = _.polyCircumradius(this.opts.bondLength, n.getSize()),
          i = r * r;
        if (e.distanceSq(n.center) < i) return !0;
      }
      return !1;
    }
    isEdgeInRing(e) {
      let t = this.graph.vertices[e.sourceId],
        n = this.graph.vertices[e.targetId];
      return this.areVerticesInSameRing(t, n);
    }
    isEdgeRotatable(e) {
      let t = this.graph.vertices[e.sourceId],
        n = this.graph.vertices[e.targetId];
      return !(
        e.bondType !== `-` ||
        t.isTerminal() ||
        n.isTerminal() ||
        (t.value.rings.length > 0 &&
          n.value.rings.length > 0 &&
          this.areVerticesInSameRing(t, n))
      );
    }
    isRingAromatic(e) {
      for (let t = 0; t < e.members.length; t++)
        if (!this.graph.vertices[e.members[t]].value.isPartOfAromaticRing)
          return !1;
      return !0;
    }
    isRingRegularPolygon(e, t = 1.15) {
      if (e.members.length < 3) return !1;
      let n = e.getPolygon(this.graph.vertices),
        r = e.center,
        i = 1 / 0,
        a = -1 / 0;
      for (let e = 0; e < n.length; e++) {
        let t = n[e].distance(r);
        (t < i && (i = t), t > a && (a = t));
      }
      return i < 1e-6 ? !1 : a / i < t;
    }
    getEdgeNormals(e) {
      let t = this.graph.vertices[e.sourceId].position,
        n = this.graph.vertices[e.targetId].position;
      return h.units(t, n);
    }
    getNonRingNeighbours(e) {
      let t = [],
        n = this.graph.vertices[e],
        r = n.neighbours;
      for (let e = 0; e < r.length; e++) {
        let i = this.graph.vertices[r[e]];
        f.intersection(n.value.rings, i.value.rings).length === 0 &&
          i.value.isBridge == 0 &&
          t.push(i);
      }
      return t;
    }
    getMinimumNonBondedDistance() {
      let e = 1 / 0;
      for (let t = 0; t < this.graph.vertices.length; t++) {
        let n = this.graph.vertices[t];
        if (n.value.isDrawn)
          for (let r = t + 1; r < this.graph.vertices.length; r++) {
            if (this.graph.hasEdge(t, r)) continue;
            let i = this.graph.vertices[r];
            if (!i.value.isDrawn) continue;
            let a = n.position.distance(i.position);
            a < e && (e = a);
          }
      }
      return e;
    }
    getRingExternalConnectionCount(e) {
      let t = this.getRing(e);
      if (!t) return 0;
      let n = new Set(t.members),
        r = new Set();
      for (let e = 0; e < t.members.length; e++) {
        let i = t.members[e],
          a = this.graph.vertices[i];
        for (let e = 0; e < a.neighbours.length; e++) {
          let t = a.neighbours[e];
          n.has(t) || r.add(`${i}:${t}`);
        }
      }
      return r.size;
    }
    resolveRigidRingOverlaps() {
      let e = this.getOverlapScore().total,
        t = this.getMinimumNonBondedDistance(),
        n = this.opts.bondLength * 0.3;
      for (let r = 0; r < this.graph.edges.length; r++) {
        let i = this.graph.edges[r];
        if (!this.isEdgeRotatable(i)) continue;
        let a = this.graph.getTreeDepth(i.sourceId, i.targetId),
          o = this.graph.getTreeDepth(i.targetId, i.sourceId),
          s = i.targetId,
          c = i.sourceId;
        a > o && ((s = i.sourceId), (c = i.targetId));
        let l = this.graph.vertices[s],
          u = this.graph.vertices[c],
          d = u.getNeighbours(s);
        if (d.length !== 2) continue;
        let f = this.graph.vertices[d[0]],
          p = this.graph.vertices[d[1]];
        if (
          f.value.rings.length !== 1 ||
          p.value.rings.length !== 1 ||
          f.value.rings[0] !== p.value.rings[0]
        )
          continue;
        let m = this.getRing(f.value.rings[0]);
        if (!m) continue;
        let h = 0,
          g = e,
          v = t,
          y = _.centralAngle(m.getSize()),
          b = Math.max(1, Math.floor(m.getSize() / 2));
        for (let e = 1; e <= b; e++) {
          let t = y * e;
          for (let e = 0; e < 2; e++) {
            let r = e === 0 ? t : -t;
            this.rotateSubtree(u.id, l.id, r, u.position);
            let i = this.getOverlapScore().total,
              a = this.getMinimumNonBondedDistance();
            (this.rotateSubtree(u.id, l.id, -r, u.position),
              !(a <= n) &&
                (i < g - 1e-6 || (Math.abs(i - g) <= 1e-6 && a > v + 1e-6)) &&
                ((h = r), (g = i), (v = a)));
          }
        }
        h !== 0 &&
          (this.rotateSubtree(u.id, l.id, h, u.position), (e = g), (t = v));
      }
      this.totalOverlapScore = e;
    }
    annotateStereochemistry() {
      for (let e = 0; e < this.graph.vertices.length; e++) {
        let t = this.graph.vertices[e];
        if (!t.value.isStereoCenter) continue;
        let n = t.neighbours;
        if (n.length + t.value.countImplicitHydrogens() !== 4) {
          t.value.isStereoCenter = !1;
          continue;
        }
        let r = T.getOrderArray(this.graph, t);
        if (r === void 0) {
          t.value.isStereoCenter = !1;
          continue;
        }
        let i =
          _.parityOfPermutation(r) ===
          (t.value.bracket.chirality === `@` ? -1 : 1)
            ? `R`
            : `S`;
        t.value.chirality = i;
        let a = r
            .map((e) => {
              let r = n[e],
                i = this.graph.vertices[r],
                a = 0;
              return (
                (a -= i.value.isStereoCenter ? 1e6 : 0),
                (a -= this.areVerticesInSameRing(i, t) ? 1e5 : 0),
                (a += i.value.isDrawn ? 1e4 : 0),
                (a += i.value.element === `C` ? 0 : 100),
                (a -= this.graph.getTreeDepth(r, t.id)),
                [a, r]
              );
            })
            .sort((e, t) => t[0] - e[0])[0][1],
          o = this._computeWedgeDirection(t, a, r, n, i);
        ((this.graph.getEdge(t.id, a).wedge = o),
          (this.graph.vertices[a].value.isDrawn = !0));
      }
    }
    _computeWedgeDirection(e, t, n, r, i) {
      let a = r.length,
        o = 0;
      for (let e = 0; e < a; e++)
        if (r[n[e]] === t) {
          o = e;
          break;
        }
      let s = [];
      for (let e = 0; e < a; e++)
        r[n[e]] !== t && s.push(this.graph.vertices[r[n[e]]].position);
      if (s.length === 2) {
        let n = this.graph.vertices[t].position,
          r = (n.x + s[0].x + s[1].x) / 3,
          i = (n.y + s[0].y + s[1].y) / 3;
        s.push({ x: 2 * e.position.x - r, y: 2 * e.position.y - i });
      }
      let c =
        (s[1].x - s[0].x) * (s[2].y - s[0].y) -
        (s[2].x - s[0].x) * (s[1].y - s[0].y);
      return (o % 2 == 0 ? c > 0 : c < 0) == (i === `R`) ? `up` : `down`;
    }
    initPseudoElements() {
      for (let e = 0; e < this.graph.vertices.length; e++) {
        let t = this.graph.vertices[e],
          n = t.neighbours,
          r = Array(n.length);
        for (let e = 0; e < n.length; e++) r[e] = this.graph.vertices[n[e]];
        if (
          t.getNeighbourCount() < 3 ||
          t.value.rings.length > 0 ||
          t.value.element === `P` ||
          (t.value.element === `C` &&
            r.length === 3 &&
            r[0].value.element === `N` &&
            r[1].value.element === `N` &&
            r[2].value.element === `N`)
        )
          continue;
        let i = 0,
          a = 0;
        for (let e = 0; e < r.length; e++) {
          let t = r[e],
            n = t.value.element,
            o = t.getNeighbourCount();
          (n !== `C` && n !== `H` && o === 1 && i++, o > 1 && a++);
        }
        if (a > 1 || i < 2) continue;
        let o = null;
        for (let e = 0; e < r.length; e++) {
          let t = r[e];
          t.getNeighbourCount() > 1 && (o = t);
        }
        for (let e = 0; e < r.length; e++) {
          let n = r[e];
          if (n.getNeighbourCount() > 1) continue;
          n.value.isDrawn = !1;
          let i = n.value.countImplicitHydrogens(),
            a = ``;
          (n.value.bracket && (a = n.value.bracket.charge || 0),
            t.value.attachPseudoElement(
              n.value.element,
              o ? o.value.element : null,
              i,
              a,
            ));
        }
      }
    }
  };
function se(e) {
  function t(e, t) {
    let n = e.length;
    if (n === 0 || (n > 0 && n - 1 in e))
      for (let r = 0; r < n && t.call(e[r], r, e[r]) !== !1; r++);
    else for (let n in e) if (t.call(e[n], n, e[n]) === !1) break;
  }
  function n(e) {
    let t = parseInt(e).toString(16);
    return t.length == 1 ? `0` + t : t;
  }
  function r(e, t, r, i) {
    return (
      (i = parseInt(i)),
      i === void 0 || i === 255
        ? `#` + n(e) + n(t) + n(r)
        : i !== 0 && `rgba(` + e + `,` + t + `,` + r + `,` + i / 255 + `)`
    );
  }
  function i(e, t, n) {
    return `M` + e + ` ` + t + `h` + n;
  }
  function a(e, t) {
    return (
      `<path stroke="` +
      e +
      `" d="` +
      t +
      `" />
`
    );
  }
  function o(e) {
    let n = ``;
    return (
      t(e, function (e, o) {
        if (((e = r.apply(null, e.split(`,`))), e === !1)) return;
        let s = [],
          c,
          l = 1;
        (t(o, function (e, t) {
          c && t[1] === c[1] && t[0] === c[0] + l
            ? l++
            : (c && (s.push(i(c[0], c[1], l)), (l = 1)), (c = t));
        }),
          s.push(i(c[0], c[1], l)),
          (n += a(e, s.join(``))));
      }),
      n
    );
  }
  function s(e) {
    let t = {},
      n = e.data,
      r = n.length,
      i = e.width,
      a = 0,
      o = 0,
      s;
    for (let e = 0; e < r; e += 4)
      n[e + 3] > 0 &&
        ((s = n[e] + `,` + n[e + 1] + `,` + n[e + 2] + `,` + n[e + 3]),
        (t[s] = t[s] || []),
        (a = (e / 4) % i),
        (o = Math.floor(e / 4 / i)),
        t[s].push([a, o]));
    return t;
  }
  let c = o(s(e)),
    l =
      `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 -0.5 ` +
      e.width +
      ` ` +
      e.height +
      `" shape-rendering="crispEdges"><g shape-rendering="crispEdges">` +
      c +
      `</g></svg>`,
    u = document.createElement(`div`);
  return ((u.innerHTML = l), u.firstElementChild);
}
var ce = l(d()),
  le = class {
    constructor(e, t, n, r, i = 0.3, a = 0, o = null, s = 1, c = !1) {
      ((this.points = e),
        (this.weights = t),
        (this.width = n),
        (this.height = r),
        (this.sigma = i),
        (this.interval = a),
        (this.opacity = s),
        (this.normalized = c),
        o === null &&
          (o = [
            `#c51b7d`,
            `#de77ae`,
            `#f1b6da`,
            `#fde0ef`,
            `#ffffff`,
            `#e6f5d0`,
            `#b8e186`,
            `#7fbc41`,
            `#4d9221`,
          ]),
        (this.colormap = o),
        (this.canvas = document.createElement(`canvas`)),
        (this.context = this.canvas.getContext(`2d`)),
        (this.canvas.width = this.width),
        (this.canvas.height = this.height));
    }
    setFromArray(e, t) {
      ((this.points = []),
        e.forEach((e) => {
          this.points.push(new h(e[0], e[1]));
        }),
        (this.weights = []),
        t.forEach((e) => {
          this.weights.push(e);
        }));
    }
    draw() {
      let e = [];
      for (let t = 0; t < this.width; t++) {
        let t = [];
        for (let e = 0; e < this.height; e++) t.push(0);
        e.push(t);
      }
      let t = 1 / (2 * this.sigma * this.sigma);
      for (let n = 0; n < this.points.length; n++) {
        let r = this.points[n],
          i = this.weights[n];
        for (let n = 0; n < this.width; n++)
          for (let a = 0; a < this.height; a++) {
            let o = n - r.x,
              s = a - r.y,
              c = (o * o + s * s) * t,
              l = i * Math.exp(-c);
            e[n][a] += l;
          }
      }
      let n = 1;
      if (!this.normalized) {
        let t = -(2 ** 53 - 1),
          r = 2 ** 53 - 1;
        for (let n = 0; n < this.width; n++)
          for (let i = 0; i < this.height; i++)
            (e[n][i] < r && (r = e[n][i]), e[n][i] > t && (t = e[n][i]));
        n = Math.max(Math.abs(r), Math.abs(t));
      }
      let r = ce.default.scale(this.colormap).domain([-1, 1]);
      for (let t = 0; t < this.width; t++)
        for (let i = 0; i < this.height; i++) {
          (this.normalized || (e[t][i] = e[t][i] / n),
            this.interval !== 0 &&
              (e[t][i] = Math.round(e[t][i] / this.interval) * this.interval));
          let [a, o, s] = r(e[t][i]).rgb();
          this.setPixel(new h(t, i), a, o, s);
        }
    }
    getImage(e) {
      let t = new Image();
      ((t.onload = () => {
        ((this.context.imageSmoothingEnabled = !1),
          this.context.drawImage(t, 0, 0, this.width, this.height),
          e && e(t));
      }),
        (t.onerror = (e) => {
          console.log(e);
        }),
        (t.src = this.canvas.toDataURL()));
    }
    getSVG() {
      return se(this.context.getImageData(0, 0, this.width, this.height));
    }
    setPixel(e, t, n, r) {
      ((this.context.fillStyle =
        `rgba(` + t + `,` + n + `,` + r + `,` + this.opacity + `)`),
        this.context.fillRect(e.x, e.y, 1, 1));
    }
  };
function ue(e) {
  let t = ``;
  for (let n = 0; n < e; n++)
    t +=
      `ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789`.charAt(
        Math.floor(Math.random() * 62),
      );
  return t;
}
var L = class e {
    constructor(e, t, n, r = !0) {
      if (
        ((this.svg =
          typeof t == `string` || t instanceof String
            ? document.getElementById(t)
            : t),
        (this.container = null),
        (this.opts = n),
        (this.uid = ue(5)),
        (this.gradientId = 0),
        (this.defaultGradients = new Map()),
        (this.backgroundItems = []),
        (this.paths = []),
        (this.vertices = []),
        (this.gradients = []),
        (this.highlights = []),
        (this.drawingWidth = 0),
        (this.drawingHeight = 0),
        (this.halfBondThickness = this.opts.bondThickness / 2),
        (this.themeManager = e),
        (this.maskElements = []),
        (this.maxX = -Number.MAX_VALUE),
        (this.maxY = -Number.MAX_VALUE),
        (this.minX = Number.MAX_VALUE),
        (this.minY = Number.MAX_VALUE),
        r)
      )
        for (; this.svg.firstChild;) this.svg.removeChild(this.svg.firstChild);
      ((this.style = document.createElementNS(
        `http://www.w3.org/2000/svg`,
        `style`,
      )),
        this.style.appendChild(
          document.createTextNode(`
                .element {
                    font: ${this.opts.fontSizeLarge}pt ${this.opts.fontFamily};
                }
                .sub {
                    font: ${this.opts.fontSizeSmall}pt ${this.opts.fontFamily};
                }
            `),
        ),
        this.svg
          ? this.svg.appendChild(this.style)
          : ((this.container = document.createElementNS(
              `http://www.w3.org/2000/svg`,
              `g`,
            )),
            this.container.appendChild(this.style)));
    }
    constructSvg() {
      let e = document.createElementNS(`http://www.w3.org/2000/svg`, `defs`),
        t = document.createElementNS(`http://www.w3.org/2000/svg`, `mask`),
        n = document.createElementNS(`http://www.w3.org/2000/svg`, `g`),
        r = document.createElementNS(`http://www.w3.org/2000/svg`, `g`),
        i = document.createElementNS(`http://www.w3.org/2000/svg`, `g`),
        a = document.createElementNS(`http://www.w3.org/2000/svg`, `g`),
        o = this.paths;
      {
        let e = document.createElementNS(`http://www.w3.org/2000/svg`, `rect`);
        (e.setAttributeNS(null, `x`, this.minX),
          e.setAttributeNS(null, `y`, this.minY),
          e.setAttributeNS(null, `width`, this.maxX - this.minX),
          e.setAttributeNS(null, `height`, this.maxY - this.minY),
          e.setAttributeNS(null, `fill`, `white`),
          t.appendChild(e));
      }
      (t.setAttributeNS(null, `id`, this.uid + `-text-mask`),
        t.setAttributeNS(null, `maskUnits`, `userSpaceOnUse`),
        t.setAttributeNS(null, `x`, this.minX),
        t.setAttributeNS(null, `y`, this.minY),
        t.setAttributeNS(null, `width`, this.maxX - this.minX),
        t.setAttributeNS(null, `height`, this.maxY - this.minY));
      for (let e of o) i.appendChild(e);
      for (let e of this.backgroundItems) n.appendChild(e);
      for (let e of this.highlights) r.appendChild(e);
      for (let e of this.vertices) a.appendChild(e);
      for (let e of this.maskElements) t.appendChild(e);
      for (let t of this.gradients) e.appendChild(t);
      if (
        (i.setAttributeNS(null, `mask`, `url(#` + this.uid + `-text-mask)`),
        this.updateViewbox(this.opts.scale),
        this.svg)
      )
        (this.svg.appendChild(e),
          this.svg.appendChild(t),
          this.svg.appendChild(n),
          this.svg.appendChild(r),
          this.svg.appendChild(i),
          this.svg.appendChild(a));
      else
        return (
          this.container.appendChild(e),
          this.container.appendChild(t),
          this.container.appendChild(n),
          this.container.appendChild(i),
          this.container.appendChild(a),
          this.container
        );
    }
    addLayer(e) {
      this.backgroundItems.push(e.firstChild);
    }
    getBondColor(e, t) {
      let n = this.themeManager.getColor(e.value.element),
        r = this.themeManager.getColor(t.value.element);
      if (n === r) return n;
      let i = this.themeManager.getColor(`C`),
        a = e.position.distance(t.position) / this.opts.bondLength;
      if (a > 0.9 && a < 1.1) {
        if (n === i) return this.getDefaultGradient(t, r, i);
        if (r === i) return this.getDefaultGradient(e, n, i);
      }
      return this.createLinearGradient(e.position, t.position, [
        [`20%`, n],
        [`80%`, r],
      ]);
    }
    getDefaultGradient(e, t, n) {
      let r = this.defaultGradients.get(e);
      return (
        r ||
        ((r = this.createRadialGradient(e.position, this.opts.bondLength, [
          [`20%`, t || this.themeManager.getColor(e.value.element)],
          [`80%`, n || this.themeManager.getColor(`C`)],
        ])),
        this.defaultGradients.set(e, r),
        r)
      );
    }
    createRadialGradient(e, t, n) {
      let r = document.createElementNS(
          `http://www.w3.org/2000/svg`,
          `radialGradient`,
        ),
        i = this.uid + `-atom-${this.gradientId++}`;
      (r.setAttributeNS(null, `id`, i),
        r.setAttributeNS(null, `gradientUnits`, `userSpaceOnUse`),
        r.setAttributeNS(null, `cx`, e.x),
        r.setAttributeNS(null, `cy`, e.y),
        r.setAttributeNS(null, `r`, t));
      for (let [e, t] of n) {
        let n = document.createElementNS(`http://www.w3.org/2000/svg`, `stop`);
        (n.setAttributeNS(null, `offset`, e),
          n.setAttributeNS(null, `stop-color`, t),
          r.appendChild(n));
      }
      return (this.gradients.push(r), `url(#${i})`);
    }
    createLinearGradient(e, t, n) {
      let r = document.createElementNS(
          `http://www.w3.org/2000/svg`,
          `linearGradient`,
        ),
        i = this.uid + `-bond-${this.gradientId++}`;
      (r.setAttributeNS(null, `id`, i),
        r.setAttributeNS(null, `gradientUnits`, `userSpaceOnUse`),
        r.setAttributeNS(null, `x1`, e.x),
        r.setAttributeNS(null, `y1`, e.y),
        r.setAttributeNS(null, `x2`, t.x),
        r.setAttributeNS(null, `y2`, t.y));
      for (let [e, t] of n) {
        let n = document.createElementNS(`http://www.w3.org/2000/svg`, `stop`);
        (n.setAttributeNS(null, `offset`, e),
          n.setAttributeNS(null, `stop-color`, t),
          r.appendChild(n));
      }
      return (this.gradients.push(r), `url(#${i})`);
    }
    createSubSuperScripts(e, t) {
      let n = document.createElementNS(`http://www.w3.org/2000/svg`, `tspan`);
      return (
        n.setAttributeNS(null, `baseline-shift`, t),
        n.appendChild(document.createTextNode(e)),
        n.setAttributeNS(null, `class`, `sub`),
        n
      );
    }
    static createUnicodeCharge(t) {
      return t === 1
        ? `⁺`
        : t === -1
          ? `⁻`
          : t > 1
            ? e.createUnicodeSuperscript(t) + `⁺`
            : t < -1
              ? e.createUnicodeSuperscript(t) + `⁻`
              : ``;
    }
    determineDimensions(e) {
      for (let t = 0; t < e.length; t++) {
        if (!e[t].value.isDrawn) continue;
        let n = e[t].position;
        (this.maxX < n.x && (this.maxX = n.x),
          this.maxY < n.y && (this.maxY = n.y),
          this.minX > n.x && (this.minX = n.x),
          this.minY > n.y && (this.minY = n.y));
      }
      let t = this.opts.padding;
      ((this.maxX += t),
        (this.maxY += t),
        (this.minX -= t),
        (this.minY -= t),
        (this.drawingWidth = this.maxX - this.minX),
        (this.drawingHeight = this.maxY - this.minY));
    }
    updateViewbox(e) {
      let t = this.minX,
        n = this.minY,
        r = this.maxX - this.minX,
        i = this.maxY - this.minY;
      if (e <= 0) {
        if (r > i) {
          let e = r - i;
          ((i = r), (n -= e / 2));
        } else {
          let e = i - r;
          ((r = i), (t -= e / 2));
        }
      } else
        this.svg &&
          ((this.svg.style.width = e * r + `px`),
          (this.svg.style.height = e * i + `px`));
      this.svg.setAttributeNS(null, `viewBox`, `${t} ${n} ${r} ${i}`);
    }
    drawBall(e, t, n) {
      let r = this.opts.bondLength / 4.5;
      (e - r < this.minX && (this.minX = e - r),
        e + r > this.maxX && (this.maxX = e + r),
        t - r < this.minY && (this.minY = t - r),
        t + r > this.maxY && (this.maxY = t + r));
      let i = document.createElementNS(`http://www.w3.org/2000/svg`, `circle`);
      (i.setAttributeNS(null, `cx`, e),
        i.setAttributeNS(null, `cy`, t),
        i.setAttributeNS(null, `r`, r),
        i.setAttributeNS(null, `fill`, this.themeManager.getColor(n)),
        this.vertices.push(i));
    }
    drawWedge(e, t = `#000`) {
      let n = e.getLeftVector().clone(),
        r = e.getRightVector().clone(),
        i = h.normals(n, r);
      (i[0].normalize(), i[1].normalize());
      let a = e.getRightChiral(),
        o = n,
        s = r;
      a && ((o = r), (s = n));
      let c = h.add(o, h.multiplyScalar(i[0], this.halfBondThickness)),
        l = h.add(s, h.multiplyScalar(i[0], 3 + this.opts.fontSizeLarge / 4)),
        u = h.add(s, h.multiplyScalar(i[1], 3 + this.opts.fontSizeLarge / 4)),
        d = h.add(o, h.multiplyScalar(i[1], this.halfBondThickness)),
        f = document.createElementNS(`http://www.w3.org/2000/svg`, `polygon`);
      (f.setAttributeNS(
        null,
        `points`,
        `${c.x},${c.y} ${l.x},${l.y} ${u.x},${u.y} ${d.x},${d.y}`,
      ),
        f.setAttributeNS(null, `fill`, t),
        this.paths.push(f));
    }
    drawAtomHighlight(e, t, n = `#03fc9d`) {
      let r = document.createElementNS(`http://www.w3.org/2000/svg`, `circle`);
      (r.setAttributeNS(null, `cx`, e),
        r.setAttributeNS(null, `cy`, t),
        r.setAttributeNS(null, `r`, this.opts.bondLength / 3),
        r.setAttributeNS(null, `fill`, n),
        this.highlights.push(r));
    }
    drawDashedWedge(e, t = `#000`) {
      if (
        isNaN(e.from.x) ||
        isNaN(e.from.y) ||
        isNaN(e.to.x) ||
        isNaN(e.to.y)
      ) {
        console.error(
          `Invalid line passed to SvgWrapper.drawDashedWedge()!`,
          e,
        );
        return;
      }
      let n = e.getLeftVector().clone(),
        r = e.getRightVector().clone(),
        i = h.normals(n, r);
      (i[0].normalize(), i[1].normalize());
      let a = e.getRightChiral(),
        o,
        s;
      a ? ((o = r), (s = n)) : ((o = n), (s = r));
      let c = h.subtract(s, o).normalize(),
        l = e.getLength(),
        u = 1.25 / (l / (this.opts.bondLength / 10));
      for (let e = 0; e < 1; e += u) {
        let n = h.multiplyScalar(c, e * l),
          r = h.add(o, n),
          a = (this.opts.fontSizeLarge / 2) * e,
          s = h.multiplyScalar(i[0], a);
        r.subtract(s);
        let u = r.clone();
        (u.add(h.multiplyScalar(s, 2)), this.drawLine(new g(r, u), !1, t));
      }
    }
    drawDebugPoint(e, t, n = ``, r = `#f00`) {
      let i = document.createElementNS(`http://www.w3.org/2000/svg`, `circle`);
      (i.setAttributeNS(null, `cx`, e),
        i.setAttributeNS(null, `cy`, t),
        i.setAttributeNS(null, `r`, `2`),
        i.setAttributeNS(null, `fill`, r),
        this.vertices.push(i),
        this.drawDebugText(e + 2, t - 2, n, r));
    }
    drawDebugText(e, t, n, r = `#f00`) {
      let i = document.createElementNS(`http://www.w3.org/2000/svg`, `text`);
      (i.setAttributeNS(null, `x`, e),
        i.setAttributeNS(null, `y`, t),
        i.setAttributeNS(null, `class`, `debug`),
        i.setAttributeNS(null, `fill`, r),
        i.setAttributeNS(null, `style`, `font: 5px sans-serif`),
        i.appendChild(document.createTextNode(n)),
        this.vertices.push(i));
    }
    drawRing(e, t, n) {
      let r = document.createElementNS(`http://www.w3.org/2000/svg`, `circle`),
        i = _.apothemFromSideLength(this.opts.bondLength, n);
      (r.setAttributeNS(null, `cx`, e),
        r.setAttributeNS(null, `cy`, t),
        r.setAttributeNS(null, `r`, i - this.opts.bondSpacing),
        r.setAttributeNS(null, `stroke`, this.themeManager.getColor(`C`)),
        r.setAttributeNS(null, `stroke-width`, this.opts.bondThickness),
        r.setAttributeNS(null, `fill`, `none`),
        this.paths.push(r));
    }
    drawLine(e, t = !1, n = `#000`, r = `round`) {
      let i = e.getLeftVector(),
        a = e.getRightVector(),
        o = document.createElementNS(`http://www.w3.org/2000/svg`, `line`);
      (o.setAttributeNS(null, `x1`, i.x),
        o.setAttributeNS(null, `y1`, i.y),
        o.setAttributeNS(null, `x2`, a.x),
        o.setAttributeNS(null, `y2`, a.y),
        o.setAttributeNS(null, `stroke`, n),
        o.setAttributeNS(null, `stroke-width`, this.opts.bondThickness),
        o.setAttributeNS(null, `stroke-linecap`, r),
        t && o.setAttributeNS(null, `stroke-dasharray`, `5,5`),
        this.paths.push(o));
    }
    drawPoint(e, t, n) {
      let r = 0.75;
      (e - r < this.minX && (this.minX = e - r),
        e + r > this.maxX && (this.maxX = e + r),
        t - r < this.minY && (this.minY = t - r),
        t + r > this.maxY && (this.maxY = t + r));
      let i = document.createElementNS(`http://www.w3.org/2000/svg`, `circle`);
      (i.setAttributeNS(null, `cx`, e),
        i.setAttributeNS(null, `cy`, t),
        i.setAttributeNS(null, `r`, `1.5`),
        i.setAttributeNS(null, `fill`, `black`),
        this.maskElements.push(i));
      let a = document.createElementNS(`http://www.w3.org/2000/svg`, `circle`);
      (a.setAttributeNS(null, `cx`, e),
        a.setAttributeNS(null, `cy`, t),
        a.setAttributeNS(null, `r`, r),
        a.setAttributeNS(null, `fill`, this.themeManager.getColor(n)),
        this.vertices.push(a));
    }
    drawText(t, n, r, i, a, o, s, c, l, u = {}) {
      let d = [],
        f = r;
      (s !== 0 && s !== null && (f += e.createUnicodeCharge(s)),
        c !== 0 && c !== null && (f = e.createUnicodeSuperscript(c) + f),
        d.push([f, r]),
        i === 1
          ? d.push([`H`, `H`])
          : i > 1 && d.push([`H` + e.createUnicodeSubscript(i), `H`]),
        s === 1 &&
          r === `N` &&
          `0O` in u &&
          `0O-1` in u &&
          ((u = {
            "0O": {
              element: `O`,
              count: 2,
              hydrogenCount: 0,
              previousElement: `C`,
              charge: ``,
            },
          }),
          (s = 0)));
      for (let t of Object.keys(u)) {
        let n = u[t],
          r = n.element;
        (n.count > 1 && (r += e.createUnicodeSubscript(n.count)),
          n.charge && (r += e.createUnicodeCharge(n.charge)),
          d.push([r, n.element]));
        let i = n.hydrogenCount * n.count;
        i === 1
          ? d.push([`H`, `H`])
          : i > 1 && d.push([`H` + e.createUnicodeSubscript(i), `H`]);
      }
      this.write(d, a, t, n, l === 1);
    }
    drawWeights(e, t, n) {
      let r = this.opts.weights.sigma,
        i = this.opts.weights.colormap,
        a = i == null ? `#4d9221` : i[i.length - 1],
        o = i == null ? `#c51b7d` : i[0],
        s = document.createElementNS(`http://www.w3.org/2000/svg`, `g`);
      (s.setAttributeNS(null, `style`, `filter:blur(${r / 2}px)`),
        this.backgroundItems.push(s));
      for (let i = 0; i < e.length; ++i) {
        let c = e[i],
          l = t[i];
        if (!c) continue;
        let u = c > 0 ? a : o;
        u = `rgb(from ${u} r g b / ${Math.abs(c) * n})`;
        let d = document.createElementNS(
          `http://www.w3.org/2000/svg`,
          `circle`,
        );
        (d.setAttributeNS(null, `cx`, l.x),
          d.setAttributeNS(null, `cy`, l.y),
          d.setAttributeNS(null, `r`, r),
          d.setAttributeNS(null, `fill`, u),
          s.appendChild(d));
      }
    }
    write(t, n, r, i, a) {
      let o = e.measureText(
        t[0][1],
        this.opts.fontSizeLarge,
        this.opts.fontFamily,
      );
      (n === `left` &&
        t[0][0] !== t[0][1] &&
        (o.width = e.measureText(
          t[0][0],
          this.opts.fontSizeLarge,
          this.opts.fontFamily,
        ).width),
        a
          ? (r + o.width * t.length > this.maxX &&
              (this.maxX = r + o.width * t.length),
            r - o.width / 2 < this.minX && (this.minX = r - o.width / 2),
            i - o.height < this.minY && (this.minY = i - o.height),
            i + o.height > this.maxY && (this.maxY = i + o.height))
          : (n === `right`
              ? n !== `left` &&
                (r + o.width * t.length > this.maxX &&
                  (this.maxX = r + o.width * t.length),
                r - o.width / 2 < this.minX && (this.minX = r - o.width / 2))
              : (r + o.width * t.length > this.maxX &&
                  (this.maxX = r + o.width * t.length),
                r - o.width * t.length < this.minX &&
                  (this.minX = r - o.width * t.length)),
            i - o.height < this.minY && (this.minY = i - o.height),
            i + o.height > this.maxY && (this.maxY = i + o.height),
            n === `down` &&
              i + 0.8 * o.height * t.length > this.maxY &&
              (this.maxY = i + 0.8 * o.height * t.length),
            n === `up` &&
              i - 0.8 * o.height * t.length < this.minY &&
              (this.minY = i - 0.8 * o.height * t.length)));
      let s = r,
        c = i,
        l = document.createElementNS(`http://www.w3.org/2000/svg`, `text`);
      l.setAttributeNS(null, `class`, `element`);
      let u = document.createElementNS(`http://www.w3.org/2000/svg`, `g`);
      (l.setAttributeNS(null, `fill`, `#ffffff`),
        n === `left` && (t = t.reverse()),
        (n === `right` || n === `down` || n === `up`) && (r -= o.width / 2),
        n === `left` && (r += o.width / 2),
        t.forEach((e, t) => {
          let r = e[0],
            i = e[1],
            a = document.createElementNS(`http://www.w3.org/2000/svg`, `tspan`);
          (a.setAttributeNS(null, `fill`, this.themeManager.getColor(i)),
            (a.textContent = r),
            (n === `up` || n === `down`) &&
              (a.setAttributeNS(null, `x`, `0px`),
              n === `up`
                ? a.setAttributeNS(null, `y`, `-${0.9 * t}em`)
                : a.setAttributeNS(null, `y`, `${0.9 * t}em`)),
            l.appendChild(a));
        }),
        l.setAttributeNS(null, `data-direction`, n),
        n === `left` || n === `right`
          ? (l.setAttributeNS(null, `dominant-baseline`, `alphabetic`),
            l.setAttributeNS(null, `y`, `0.36em`))
          : l.setAttributeNS(null, `dominant-baseline`, `central`),
        n === `left` && l.setAttributeNS(null, `text-anchor`, `end`),
        u.appendChild(l),
        u.setAttributeNS(
          null,
          `style`,
          `transform: translateX(${r}px) translateY(${i}px)`,
        ));
      let d = this.opts.fontSizeLarge * 0.75;
      t[0][1].length > 1 && (d = this.opts.fontSizeLarge * 1.1);
      let f = document.createElementNS(`http://www.w3.org/2000/svg`, `circle`);
      (f.setAttributeNS(null, `cx`, s),
        f.setAttributeNS(null, `cy`, c),
        f.setAttributeNS(null, `r`, d),
        f.setAttributeNS(null, `fill`, `black`),
        this.maskElements.push(f),
        this.vertices.push(u));
    }
    toCanvas(e, t, n) {
      (typeof e == `string` || e instanceof String) &&
        (e = document.getElementById(e));
      let r = new Image();
      ((r.onload = function () {
        ((e.width = t),
          (e.height = n),
          e.getContext(`2d`).drawImage(r, 0, 0, t, n));
      }),
        (r.src =
          `data:image/svg+xml;charset-utf-8,` +
          encodeURIComponent(this.svg.outerHTML)));
    }
    static createUnicodeSubscript(e) {
      let t = ``;
      return (
        e
          .toString()
          .split(``)
          .forEach((e) => {
            t += [`₀`, `₁`, `₂`, `₃`, `₄`, `₅`, `₆`, `₇`, `₈`, `₉`][
              parseInt(e)
            ];
          }),
        t
      );
    }
    static createUnicodeSuperscript(e) {
      let t = ``;
      return (
        e
          .toString()
          .split(``)
          .forEach((e) => {
            let n = parseInt(e);
            Number.isFinite(n) &&
              (t += [`⁰`, `¹`, `²`, `³`, `⁴`, `⁵`, `⁶`, `⁷`, `⁸`, `⁹`][n]);
          }),
        t
      );
    }
    static replaceNumbersWithSubscript(e) {
      for (let [t, n] of Object.entries({
        0: `₀`,
        1: `₁`,
        2: `₂`,
        3: `₃`,
        4: `₄`,
        5: `₅`,
        6: `₆`,
        7: `₇`,
        8: `₈`,
        9: `₉`,
      }))
        e = e.replaceAll(t, n);
      return e;
    }
    static measureText(t, n, r) {
      let i = document.createElement(`canvas`).getContext(`2d`);
      if (!i) return e.estimateTextSize(t, n);
      let a = n / 16;
      i.font = `16pt ${r}`;
      let o = i.measureText(t),
        s =
          Math.abs(o.actualBoundingBoxLeft) +
          Math.abs(o.actualBoundingBoxRight),
        c =
          Math.abs(o.actualBoundingBoxAscent) +
          Math.abs(o.actualBoundingBoxDescent);
      return { width: s * a, height: c * a };
    }
    static estimateTextSize(e, t) {
      let n = 0;
      for (let t of String(e))
        t === ` `
          ? (n += 0.35)
          : /[A-Z]/.test(t)
            ? (n += 0.68)
            : /[a-z]/.test(t)
              ? (n += 0.56)
              : /[0-9]/.test(t)
                ? (n += 0.55)
                : (n += 0.45);
      return { width: t * n, height: t };
    }
    static svgToCanvas(e, t, n, r, i = null) {
      (e.setAttributeNS(null, `width`, n), e.setAttributeNS(null, `height`, r));
      let a = document.createElement(`img`);
      return (
        (a.onload = function () {
          ((t.width = n), (t.height = r));
          let e = t.getContext(`2d`);
          ((e.imageSmoothingEnabled = !1),
            e.drawImage(a, 0, 0, n, r),
            i && i(t));
        }),
        (a.onerror = function (e) {
          console.log(e);
        }),
        (a.src =
          `data:image/svg+xml;charset-utf-8,` +
          encodeURIComponent(e.outerHTML)),
        t
      );
    }
    static svgToImg(e, t, n, r) {
      let i = document.createElement(`canvas`);
      this.svgToCanvas(e, i, n, r, () => {
        t.src = i.toDataURL(`image/png`);
      });
    }
    static writeText(t, n, r, i, a = 2 ** 53 - 1) {
      let o = document.createElementNS(`http://www.w3.org/2000/svg`, `svg`),
        s = document.createElementNS(`http://www.w3.org/2000/svg`, `style`);
      (s.appendChild(
        document.createTextNode(`
            .text {
                font: ${r}pt ${i};
                dominant-baseline: alphabetic;
            }
        `),
      ),
        o.appendChild(s));
      let c = document.createElementNS(`http://www.w3.org/2000/svg`, `text`);
      c.setAttributeNS(null, `class`, `text`);
      let l = [];
      t.split(
        `
`,
      ).forEach((t) => {
        if (e.measureText(t, r, i).width >= a) {
          let n = t.split(` `),
            o = 0;
          for (let t = 0; t < n.length; t++) {
            let s = n.slice(o, t + 1).join(` `);
            e.measureText(s, r, i).width > a &&
              t > o &&
              (l.push(n.slice(o, t).join(` `)), (o = t));
          }
          o < n.length && l.push(n.slice(o, n.length).join(` `));
        } else l.push(t);
      });
      let u = 0;
      return (
        l.forEach((t, a) => {
          let o = document.createElementNS(
            `http://www.w3.org/2000/svg`,
            `tspan`,
          );
          (o.setAttributeNS(null, `fill`, n.getColor(`C`)),
            (o.textContent = t),
            o.setAttributeNS(null, `x`, `0px`),
            o.setAttributeNS(null, `y`, `${a + 1}em`),
            c.appendChild(o));
          let s = e.measureText(t, r, i);
          u = Math.max(u, s.width);
        }),
        c.setAttributeNS(null, `transform`, `translate(${u / 2}, 0)`),
        c.setAttributeNS(null, `text-anchor`, `middle`),
        o.appendChild(c),
        { svg: o, width: u, height: ((l.length + 0.4) * r * 4) / 3 }
      );
    }
  },
  R = class {
    constructor(e, t = !0) {
      ((this.preprocessor = new oe(e)),
        (this.opts = this.preprocessor.opts),
        (this.clear = t),
        (this.svgWrapper = null));
    }
    draw(e, t, n = `light`, r = null, i = !1, a = [], o = !1) {
      let s = null;
      if (
        (t === null || t === `svg`
          ? ((s = document.createElementNS(
              `http://www.w3.org/2000/svg`,
              `svg`,
            )),
            s.setAttribute(`xmlns`, `http://www.w3.org/2000/svg`),
            s.setAttribute(`xmlns:xlink`, `http://www.w3.org/1999/xlink`),
            s.setAttributeNS(null, `width`, this.opts.width),
            s.setAttributeNS(null, `height`, this.opts.height))
          : (s =
              t instanceof String
                ? document.getElementById(t.valueOf())
                : typeof t == `string`
                  ? document.getElementById(t)
                  : t),
        !(s instanceof SVGSVGElement))
      )
        throw Error(`Second argument was not an SVG or the ID of an SVG.`);
      let c = {
        padding: this.opts.padding,
        compactDrawing: this.opts.compactDrawing,
      };
      r !== null &&
        ((this.opts.padding += this.opts.weights.additionalPadding),
        (this.opts.compactDrawing = !1));
      let l = this.preprocessor;
      return (
        l.initDraw(e, n, i, a),
        i ||
          ((this.themeManager = new x(this.opts.themes, n)),
          (this.svgWrapper === null || this.clear) &&
            (this.svgWrapper = new L(
              this.themeManager,
              s,
              this.opts,
              this.clear,
            ))),
        l.processGraph(),
        this.svgWrapper.determineDimensions(l.graph.vertices),
        this.drawAtomHighlights(l.opts.debug),
        this.drawEdges(l.opts.debug),
        this.drawVertices(l.opts.debug),
        r !== null && this.drawWeights(r, o),
        l.opts.debug &&
          console.debug(`SvgDrawer::draw()`, {
            graph: l.graph,
            rings: l.rings,
            ringConnections: l.ringConnections,
          }),
        this.svgWrapper.constructSvg(),
        r !== null &&
          ((this.opts.padding = c.padding),
          (this.opts.compactDrawing = c.padding)),
        s
      );
    }
    drawCanvas(e, t, n = `light`, r = !1) {
      let i = null;
      if (
        ((i =
          t instanceof String
            ? document.getElementById(t.valueOf())
            : typeof t == `string`
              ? document.getElementById(t)
              : t),
        !(i instanceof HTMLCanvasElement))
      )
        throw Error(`Second argument was not a canvas or the ID of a canvas.`);
      let a = document.createElementNS(`http://www.w3.org/2000/svg`, `svg`);
      return (
        a.setAttribute(`xmlns`, `http://www.w3.org/2000/svg`),
        a.setAttributeNS(null, `viewBox`, `0 0 500 500`),
        a.setAttributeNS(null, `width`, `500`),
        a.setAttributeNS(null, `height`, `500`),
        a.setAttributeNS(
          null,
          `style`,
          `visibility: hidden: position: absolute; left: -1000px`,
        ),
        document.body.appendChild(a),
        this.draw(e, a, n, r),
        this.svgWrapper.toCanvas(i, this.opts.width, this.opts.height),
        document.body.removeChild(a),
        t
      );
    }
    drawAromaticityRing(e) {
      this.svgWrapper.drawRing(e.center.x, e.center.y, e.getSize());
    }
    drawEdges(e) {
      let t = this.preprocessor,
        n = t.graph,
        r = t.rings,
        i = Array(this.preprocessor.graph.edges.length);
      (i.fill(!1),
        n.traverseBF(0, (t) => {
          let r = n.getEdges(t.id);
          for (let t = 0; t < r.length; t++) {
            let n = r[t];
            i[n] || ((i[n] = !0), this.drawEdge(n, e));
          }
        }));
      for (let e = 0; e < r.length; e++) {
        let n = r[e];
        t.isRingAromatic(n) &&
          ((n.isPartOfBridged && !t.isRingRegularPolygon(n)) ||
            this.drawAromaticityRing(n));
      }
    }
    drawEdge(e, t) {
      let n = this.preprocessor,
        r = n.opts,
        i = this.svgWrapper,
        a = n.graph.edges[e],
        o = n.graph.vertices[a.sourceId],
        s = n.graph.vertices[a.targetId],
        c = o.value.element,
        l = s.value.element;
      if (
        (!o.value.isDrawn || !s.value.isDrawn) &&
        n.opts.atomVisualization === "default"
      )
        return;
      let u = i.getBondColor(o, s),
        d = o.position,
        p = s.position,
        m = n.getEdgeNormals(a),
        _ = f.clone(m);
      (_[0].multiplyScalar(10).add(d), _[1].multiplyScalar(10).add(d));
      let v =
        a.isPartOfAromaticRing &&
        o.value.bridgedRing !== null &&
        s.value.bridgedRing !== null
          ? n.getLargestOrAromaticCommonRing(o, s)
          : null;
      if (
        a.bondType === `=` ||
        n.getRingbondType(o, s) === `=` ||
        (v && !n.isRingRegularPolygon(v))
      ) {
        let e = n.areVerticesInSameRing(o, s),
          t = n.chooseSide(o, s, _);
        if (e) {
          let e = n.getLargestOrAromaticCommonRing(o, s).center;
          (m[0].multiplyScalar(r.bondSpacing),
            m[1].multiplyScalar(r.bondSpacing));
          let t = null;
          ((t = e.sameSideAs(o.position, s.position, h.add(d, m[0]))
            ? new g(h.add(d, m[0]), h.add(p, m[0]), c, l)
            : new g(h.add(d, m[1]), h.add(p, m[1]), c, l)),
            t.shorten(r.bondLength - r.shortBondLength * r.bondLength),
            i.drawLine(t, a.isPartOfAromaticRing, u),
            i.drawLine(new g(d, p, c, l), !1, u));
        } else if (
          a.center ||
          (o.isTerminal() && s.isTerminal()) ||
          (t.anCount == 0 && t.bnCount > 1) ||
          (t.bnCount == 0 && t.anCount > 1)
        ) {
          this.multiplyNormals(m, r.halfBondSpacing);
          let e = new g(h.add(d, m[0]), h.add(p, m[0]), c, l),
            t = new g(h.add(d, m[1]), h.add(p, m[1]), c, l);
          (i.drawLine(e, !1, u), i.drawLine(t, !1, u));
        } else if (
          t.sideCount[0] > t.sideCount[1] ||
          t.totalSideCount[0] > t.totalSideCount[1]
        ) {
          this.multiplyNormals(m, r.bondSpacing);
          let e = new g(h.add(d, m[0]), h.add(p, m[0]), c, l);
          (e.shorten(r.bondLength - r.shortBondLength * r.bondLength),
            i.drawLine(e, !1, u),
            i.drawLine(new g(d, p, c, l), !1, u));
        } else if (
          t.sideCount[0] < t.sideCount[1] ||
          t.totalSideCount[0] <= t.totalSideCount[1]
        ) {
          this.multiplyNormals(m, r.bondSpacing);
          let e = new g(h.add(d, m[1]), h.add(p, m[1]), c, l);
          (e.shorten(r.bondLength - r.shortBondLength * r.bondLength),
            i.drawLine(e, !1, u),
            i.drawLine(new g(d, p, c, l), !1, u));
        }
      } else if (a.bondType === `#`) {
        (m[0].multiplyScalar(r.bondSpacing / 1.5),
          m[1].multiplyScalar(r.bondSpacing / 1.5));
        let e = new g(h.add(d, m[0]), h.add(p, m[0]), c, l),
          t = new g(h.add(d, m[1]), h.add(p, m[1]), c, l);
        (i.drawLine(e, !1, u),
          i.drawLine(t, !1, u),
          i.drawLine(new g(d, p, c, l), !1, u));
      } else if (a.bondType !== `.`) {
        let e = o.value.isStereoCenter,
          t = s.value.isStereoCenter;
        a.wedge === `up`
          ? i.drawWedge(new g(d, p, c, l, e, t), u)
          : a.wedge === `down`
            ? i.drawDashedWedge(new g(d, p, c, l, e, t), u)
            : i.drawLine(new g(d, p, c, l, e, t), !1, u);
      }
      if (t) {
        let t = h.midpoint(d, p);
        i.drawDebugText(t.x, t.y, `e` + e, `#0c0`);
      }
    }
    drawAtomHighlights(e) {
      let t = this.preprocessor,
        n = t.graph,
        r = this.svgWrapper;
      for (let e = 0; e < n.vertices.length; e++) {
        let i = n.vertices[e],
          a = i.value;
        for (let e = 0; e < t.highlight_atoms.length; e++) {
          let n = t.highlight_atoms[e];
          a.class === n[0] &&
            r.drawAtomHighlight(i.position.x, i.position.y, n[1]);
        }
      }
    }
    drawVertices(e) {
      let t = this.preprocessor,
        n = t.opts,
        r = t.graph,
        i = t.rings,
        a = this.svgWrapper;
      for (let t = 0; t < r.vertices.length; t++) {
        let i = r.vertices[t],
          o = i.value,
          s = 0,
          c = 0,
          l = o.element,
          u = o.countImplicitHydrogens(),
          d = i.getTextDirection(r.vertices, o.hasAttachedPseudoElements),
          p = oe.getEffectiveShowCarbonsMode(n),
          m =
            p === `terminal` || l !== `C` || o.hasAttachedPseudoElements
              ? i.isTerminal()
              : !1,
          g = o.element === `C`;
        if (l === `C`) {
          let e = o.rings && o.rings.length > 0;
          p === `none`
            ? ((g = !0), (m = !1))
            : (p === `all` || (p === `acyclic` && !e)) && ((g = !1), (m = !0));
        }
        if (
          (o.bracket && ((s = o.bracket.charge), (c = o.bracket.isotope)),
          (s || c || r.vertices.length < 3) && (g = !1),
          n.atomVisualization === `allballs`)
        )
          a.drawBall(i.position.x, i.position.y, l);
        else if (
          (o.isDrawn &&
            (!g || o.drawExplicit || m || o.hasAttachedPseudoElements)) ||
          r.vertices.length === 1
        ) {
          if (n.atomVisualization === "default") {
            let e = o.getAttachedPseudoElements();
            (o.hasAttachedPseudoElements &&
              r.vertices.length === Object.keys(e).length + 1 &&
              (d = `right`),
              a.drawText(
                i.position.x,
                i.position.y,
                l,
                u,
                d,
                m,
                s,
                c,
                r.vertices.length,
                e,
              ));
          } else
            n.atomVisualization === `balls` &&
              a.drawBall(i.position.x, i.position.y, l);
        } else if (i.getNeighbourCount() === 2 && i.forcePositioned == 1) {
          let e = r.vertices[i.neighbours[0]].position,
            t = r.vertices[i.neighbours[1]].position,
            n = h.threePointangle(i.position, e, t);
          Math.abs(Math.PI - n) < 0.1 &&
            a.drawPoint(i.position.x, i.position.y, l);
        }
        if (e) {
          let e = `v` + i.id + ` ` + f.print(o.ringbonds);
          a.drawDebugText(i.position.x, i.position.y, e);
        }
      }
      if (n.debug)
        for (let e = 0; e < i.length; e++) {
          let t = i[e].center;
          a.drawDebugPoint(t.x, t.y, `r` + i[e].id, `#00f`);
        }
    }
    drawWeights(e, t) {
      if (!e) return;
      let n = this.preprocessor.graph.atomIdxToVertexId;
      e.length < n.length
        ? (n = n.slice(0, e.length))
        : e.length > n.length &&
          (console.warn(
            `More weights (${e.length}) than heavy atoms (${n.length}); truncating.`,
          ),
          (e = e.slice(0, n.length)));
      let r = 0,
        i = 0;
      for (let t = 0; t < e.length; ++t) {
        let n = e[t];
        n && (n < r && (r = n), n > i && (i = n));
      }
      if (r === 0 && i === 0) return;
      if (this.opts.experimentalWeights) {
        let a = n.map((e) => this.preprocessor.graph.vertices[e].position),
          o = this.opts.weights.opacity;
        return (
          t || (o /= Math.max(-r, i)),
          this.svgWrapper.drawWeights(e, a, o)
        );
      }
      let a = this.svgWrapper.minX,
        o = this.svgWrapper.minY,
        s = new le(
          n.map((e) => {
            let t = this.preprocessor.graph.vertices[e];
            return new h(t.position.x - a, t.position.y - o);
          }),
          e,
          this.svgWrapper.drawingWidth,
          this.svgWrapper.drawingHeight,
          this.opts.weights.sigma,
          this.opts.weights.interval,
          this.opts.weights.colormap,
          this.opts.weights.opacity,
          t,
        );
      s.draw();
      let c = s.getSVG();
      (c.firstElementChild.setAttributeNS(
        null,
        `transform`,
        `translate(${a},${o})`,
      ),
        this.svgWrapper.addLayer(c));
    }
    getTotalOverlapScore() {
      return this.preprocessor.getTotalOverlapScore();
    }
    getMolecularFormula(e = null) {
      return this.preprocessor.getMolecularFormula(e);
    }
    multiplyNormals(e, t) {
      (e[0].multiplyScalar(t), e[1].multiplyScalar(t));
    }
  },
  z = class {
    constructor(e) {
      this.svgDrawer = new R(e);
    }
    draw(e, t, n = `light`, r = !1, i = []) {
      let a = null,
        o = null;
      if (
        ((a =
          t instanceof String
            ? document.getElementById(t.valueOf())
            : typeof t == `string`
              ? document.getElementById(t)
              : t),
        a instanceof HTMLCanvasElement)
      )
        o = a;
      else
        throw Error(`Second argument was not a canvas or the ID of a canvas.`);
      let s = document.createElementNS(`http://www.w3.org/2000/svg`, `svg`);
      (s.setAttribute(`xmlns`, `http://www.w3.org/2000/svg`),
        s.setAttributeNS(
          null,
          `viewBox`,
          `0 0 ` + this.svgDrawer.opts.width + ` ` + this.svgDrawer.opts.height,
        ),
        s.setAttributeNS(null, `width`, this.svgDrawer.opts.width + ``),
        s.setAttributeNS(null, `height`, this.svgDrawer.opts.height + ``),
        this.svgDrawer.draw(e, s, n, null, r, i),
        this.svgDrawer.svgWrapper.toCanvas(
          o,
          this.svgDrawer.opts.width,
          this.svgDrawer.opts.height,
        ));
    }
    getTotalOverlapScore() {
      return this.svgDrawer.getTotalOverlapScore();
    }
    getMolecularFormula() {
      return this.svgDrawer.getMolecularFormula();
    }
  },
  B = (function () {
    function e(e, t) {
      function n() {
        this.constructor = e;
      }
      ((n.prototype = t.prototype), (e.prototype = new n()));
    }
    function t(e, n, r, i) {
      ((this.message = e),
        (this.expected = n),
        (this.found = r),
        (this.location = i),
        (this.name = `SyntaxError`),
        typeof Error.captureStackTrace == `function` &&
          Error.captureStackTrace(this, t));
    }
    (e(t, Error),
      (t.buildMessage = function (e, t) {
        var n = {
          literal: function (e) {
            return `"` + i(e.text) + `"`;
          },
          class: function (e) {
            var t = ``,
              n;
            for (n = 0; n < e.parts.length; n++)
              t +=
                e.parts[n] instanceof Array
                  ? a(e.parts[n][0]) + `-` + a(e.parts[n][1])
                  : a(e.parts[n]);
            return `[` + (e.inverted ? `^` : ``) + t + `]`;
          },
          any: function (e) {
            return `any character`;
          },
          end: function (e) {
            return `end of input`;
          },
          other: function (e) {
            return e.description;
          },
        };
        function r(e) {
          return e.charCodeAt(0).toString(16).toUpperCase();
        }
        function i(e) {
          return e
            .replace(/\\/g, `\\\\`)
            .replace(/"/g, `\\"`)
            .replace(/\0/g, `\\0`)
            .replace(/\t/g, `\\t`)
            .replace(/\n/g, `\\n`)
            .replace(/\r/g, `\\r`)
            .replace(/[\x00-\x0F]/g, function (e) {
              return `\\x0` + r(e);
            })
            .replace(/[\x10-\x1F\x7F-\x9F]/g, function (e) {
              return `\\x` + r(e);
            });
        }
        function a(e) {
          return e
            .replace(/\\/g, `\\\\`)
            .replace(/\]/g, `\\]`)
            .replace(/\^/g, `\\^`)
            .replace(/-/g, `\\-`)
            .replace(/\0/g, `\\0`)
            .replace(/\t/g, `\\t`)
            .replace(/\n/g, `\\n`)
            .replace(/\r/g, `\\r`)
            .replace(/[\x00-\x0F]/g, function (e) {
              return `\\x0` + r(e);
            })
            .replace(/[\x10-\x1F\x7F-\x9F]/g, function (e) {
              return `\\x` + r(e);
            });
        }
        function o(e) {
          return n[e.type](e);
        }
        function s(e) {
          var t = Array(e.length),
            n,
            r;
          for (n = 0; n < e.length; n++) t[n] = o(e[n]);
          if ((t.sort(), t.length > 0)) {
            for (n = 1, r = 1; n < t.length; n++)
              t[n - 1] !== t[n] && ((t[r] = t[n]), r++);
            t.length = r;
          }
          switch (t.length) {
            case 1:
              return t[0];
            case 2:
              return t[0] + ` or ` + t[1];
            default:
              return t.slice(0, -1).join(`, `) + `, or ` + t[t.length - 1];
          }
        }
        function c(e) {
          return e ? `"` + i(e) + `"` : `end of input`;
        }
        return `Expected ` + s(e) + ` but ` + c(t) + ` found.`;
      }));
    function n(e, n) {
      if (
        ((n = n === void 0 ? {} : n),
        e.split(`(`).length - 1 != e.split(`)`).length - 1)
      )
        throw Xe(
          `The number of opening parentheses does not match the number of closing parentheses.`,
          0,
        );
      var r = {},
        i = { chain: Qe },
        a = Qe,
        o = function (e) {
          for (var t = [], n = [], r = 0; r < e[1].length; r++) t.push(e[1][r]);
          for (var r = 0; r < e[2].length; r++) {
            var i = e[2][r][0] ? e[2][r][0] : `-`;
            n.push({ bond: i, id: e[2][r][1] });
          }
          for (var r = 0; r < e[3].length; r++) t.push(e[3][r]);
          for (var r = 0; r < e[6].length; r++) t.push(e[6][r]);
          return {
            atom: e[0],
            isBracket: !!e[0].element,
            branches: t,
            branchCount: t.length,
            ringbonds: n,
            ringbondCount: n.length,
            bond: e[4] ? e[4] : `-`,
            next: e[5],
            hasNext: !!e[5],
          };
        },
        s = `(`,
        c = G(`(`, !1),
        l = `)`,
        u = G(`)`, !1),
        d = function (e) {
          var t = e[1] ? e[1] : `-`;
          return ((e[2].branchBond = t), e[2]);
        },
        f = function (e) {
          return e;
        },
        p = /^[\-=#$:\/\\.]/,
        m = K([`-`, `=`, `#`, `$`, `:`, `/`, `\\`, `.`], !1, !1),
        h = function (e) {
          return e;
        },
        g = `[`,
        _ = G(`[`, !1),
        v = `se`,
        y = G(`se`, !1),
        b = `as`,
        x = G(`as`, !1),
        S = `]`,
        C = G(`]`, !1),
        w = function (e) {
          return {
            isotope: e[1],
            element: e[2],
            chirality: e[3],
            hcount: e[4],
            charge: e[5],
            class: e[6],
          };
        },
        T = `B`,
        E = G(`B`, !1),
        D = `r`,
        O = G(`r`, !1),
        k = `C`,
        A = G(`C`, !1),
        ee = `l`,
        j = G(`l`, !1),
        te = /^[NOPSFI]/,
        M = K([`N`, `O`, `P`, `S`, `F`, `I`], !1, !1),
        N = function (e) {
          return e.length > 1 ? e.join(``) : e;
        },
        P = /^[bcnops]/,
        ne = K([`b`, `c`, `n`, `o`, `p`, `s`], !1, !1),
        F = `*`,
        I = G(`*`, !1),
        re = function (e) {
          return e;
        },
        ie = /^[A-Z]/,
        ae = K([[`A`, `Z`]], !1, !1),
        oe = /^[a-z]/,
        se = K([[`a`, `z`]], !1, !1),
        ce = function (e) {
          return e.join(``);
        },
        le = `%`,
        ue = G(`%`, !1),
        L = /^[1-9]/,
        R = K([[`1`, `9`]], !1, !1),
        z = /^[0-9]/,
        B = K([[`0`, `9`]], !1, !1),
        de = function (e) {
          return e.length == 1
            ? Number(e)
            : Number(e.join(``).replace(`%`, ``));
        },
        fe = `@`,
        pe = G(`@`, !1),
        me = `TH`,
        he = G(`TH`, !1),
        ge = /^[12]/,
        V = K([`1`, `2`], !1, !1),
        _e = `AL`,
        ve = G(`AL`, !1),
        ye = `SP`,
        be = G(`SP`, !1),
        xe = /^[1-3]/,
        Se = K([[`1`, `3`]], !1, !1),
        Ce = `TB`,
        we = G(`TB`, !1),
        Te = `OH`,
        Ee = G(`OH`, !1),
        De = function (e) {
          return e[1]
            ? e[1] == `@`
              ? `@@`
              : e[1].join(``).replace(`,`, ``)
            : `@`;
        },
        Oe = function (e) {
          return e;
        },
        ke = `+`,
        Ae = G(`+`, !1),
        je = function (e) {
          return e[1] ? (e[1] == `+` ? 2 : Number(e[1].join(``))) : 1;
        },
        Me = `-`,
        Ne = G(`-`, !1),
        Pe = function (e) {
          return e[1] ? (e[1] == `-` ? -2 : -Number(e[1].join(``))) : -1;
        },
        Fe = `H`,
        Ie = G(`H`, !1),
        Le = function (e) {
          return e[1] ? Number(e[1]) : 1;
        },
        Re = `:`,
        ze = G(`:`, !1),
        Be = /^[0]/,
        Ve = K([`0`], !1, !1),
        He = function (e) {
          return Number(e[1][0] + e[1][1].join(``));
        },
        Ue = function (e) {
          return Number(e.join(``));
        },
        H = 0,
        We = [{ line: 1, column: 1 }],
        U = 0,
        Ge = [],
        W = 0,
        Ke;
      if (`startRule` in n) {
        if (!(n.startRule in i))
          throw Error(`Can't start parsing from rule "` + n.startRule + `".`);
        a = i[n.startRule];
      }
      function G(e, t) {
        return { type: `literal`, text: e, ignoreCase: t };
      }
      function K(e, t, n) {
        return { type: `class`, parts: e, inverted: t, ignoreCase: n };
      }
      function qe() {
        return { type: `end` };
      }
      function Je(t) {
        var n = We[t],
          r;
        if (n) return n;
        for (r = t - 1; !We[r];) r--;
        for (n = We[r], n = { line: n.line, column: n.column }; r < t;)
          (e.charCodeAt(r) === 10 ? (n.line++, (n.column = 1)) : n.column++,
            r++);
        return ((We[t] = n), n);
      }
      function Ye(e, t) {
        var n = Je(e),
          r = Je(t);
        return {
          start: { offset: e, line: n.line, column: n.column },
          end: { offset: t, line: r.line, column: r.column },
        };
      }
      function q(e) {
        H < U || (H > U && ((U = H), (Ge = [])), Ge.push(e));
      }
      function Xe(e, n) {
        return new t(e, null, null, n);
      }
      function Ze(e, n, r) {
        return new t(t.buildMessage(e, n), e, n, r);
      }
      function Qe() {
        var e, t, n, i, a, s, c, l, u, d;
        if (((e = H), (t = H), (n = $e()), n !== r)) {
          for (i = [], a = J(); a !== r;) (i.push(a), (a = J()));
          if (i !== r) {
            for (
              a = [],
                s = H,
                c = et(),
                c === r && (c = null),
                c === r
                  ? ((H = s), (s = r))
                  : ((l = ot()),
                    l === r ? ((H = s), (s = r)) : ((c = [c, l]), (s = c)));
              s !== r;
            )
              (a.push(s),
                (s = H),
                (c = et()),
                c === r && (c = null),
                c === r
                  ? ((H = s), (s = r))
                  : ((l = ot()),
                    l === r ? ((H = s), (s = r)) : ((c = [c, l]), (s = c))));
            if (a !== r) {
              for (s = [], c = J(); c !== r;) (s.push(c), (c = J()));
              if (s !== r) {
                if (((c = et()), c === r && (c = null), c !== r)) {
                  if (((l = Qe()), l === r && (l = null), l !== r)) {
                    for (u = [], d = J(); d !== r;) (u.push(d), (d = J()));
                    u === r
                      ? ((H = t), (t = r))
                      : ((n = [n, i, a, s, c, l, u]), (t = n));
                  } else ((H = t), (t = r));
                } else ((H = t), (t = r));
              } else ((H = t), (t = r));
            } else ((H = t), (t = r));
          } else ((H = t), (t = r));
        } else ((H = t), (t = r));
        return (t !== r && (t = o(t)), (e = t), e);
      }
      function J() {
        var t, n, i, a, o, f;
        return (
          (t = H),
          (n = H),
          e.charCodeAt(H) === 40 ? ((i = s), H++) : ((i = r), W === 0 && q(c)),
          i === r
            ? ((H = n), (n = r))
            : ((a = et()),
              a === r && (a = null),
              a === r
                ? ((H = n), (n = r))
                : ((o = Qe()),
                  o === r
                    ? ((H = n), (n = r))
                    : (e.charCodeAt(H) === 41
                        ? ((f = l), H++)
                        : ((f = r), W === 0 && q(u)),
                      f === r
                        ? ((H = n), (n = r))
                        : ((i = [i, a, o, f]), (n = i))))),
          n !== r && (n = d(n)),
          (t = n),
          t
        );
      }
      function $e() {
        var e, t;
        return (
          (e = H),
          (t = nt()),
          t === r &&
            ((t = rt()), t === r && ((t = tt()), t === r && (t = it()))),
          t !== r && (t = f(t)),
          (e = t),
          e
        );
      }
      function et() {
        var t, n;
        if (((t = H), p.test(e.charAt(H)))) {
          if (((n = e.charAt(H)), n === e.charAt(H + 1) && ((n = r), W === 0)))
            throw Xe(`The parser encountered a bond repetition.`, H + 1);
          H++;
        } else ((n = r), W === 0 && q(m));
        return (n !== r && (n = h(n)), (t = n), t);
      }
      function tt() {
        var t, n, i, a, o, s, c, l, u, d;
        return (
          (t = H),
          (n = H),
          e.charCodeAt(H) === 91 ? ((i = g), H++) : ((i = r), W === 0 && q(_)),
          i === r
            ? ((H = n), (n = r))
            : ((a = pt()),
              a === r && (a = null),
              a === r
                ? ((H = n), (n = r))
                : (e.substr(H, 2) === v
                    ? ((o = v), (H += 2))
                    : ((o = r), W === 0 && q(y)),
                  o === r &&
                    (e.substr(H, 2) === b
                      ? ((o = b), (H += 2))
                      : ((o = r), W === 0 && q(x)),
                    o === r &&
                      ((o = rt()),
                      o === r && ((o = at()), o === r && (o = it())))),
                  o === r
                    ? ((H = n), (n = r))
                    : ((s = st()),
                      s === r && (s = null),
                      s === r
                        ? ((H = n), (n = r))
                        : ((c = dt()),
                          c === r && (c = null),
                          c === r
                            ? ((H = n), (n = r))
                            : ((l = ct()),
                              l === r && (l = null),
                              l === r
                                ? ((H = n), (n = r))
                                : ((u = ft()),
                                  u === r && (u = null),
                                  u === r
                                    ? ((H = n), (n = r))
                                    : (e.charCodeAt(H) === 93
                                        ? ((d = S), H++)
                                        : ((d = r), W === 0 && q(C)),
                                      d === r
                                        ? ((H = n), (n = r))
                                        : ((i = [i, a, o, s, c, l, u, d]),
                                          (n = i))))))))),
          n !== r && (n = w(n)),
          (t = n),
          t
        );
      }
      function nt() {
        var t, n, i, a;
        return (
          (t = H),
          (n = H),
          e.charCodeAt(H) === 66 ? ((i = T), H++) : ((i = r), W === 0 && q(E)),
          i === r
            ? ((H = n), (n = r))
            : (e.charCodeAt(H) === 114
                ? ((a = D), H++)
                : ((a = r), W === 0 && q(O)),
              a === r && (a = null),
              a === r ? ((H = n), (n = r)) : ((i = [i, a]), (n = i))),
          n === r &&
            ((n = H),
            e.charCodeAt(H) === 67
              ? ((i = k), H++)
              : ((i = r), W === 0 && q(A)),
            i === r
              ? ((H = n), (n = r))
              : (e.charCodeAt(H) === 108
                  ? ((a = ee), H++)
                  : ((a = r), W === 0 && q(j)),
                a === r && (a = null),
                a === r ? ((H = n), (n = r)) : ((i = [i, a]), (n = i))),
            n === r &&
              (te.test(e.charAt(H))
                ? ((n = e.charAt(H)), H++)
                : ((n = r), W === 0 && q(M)))),
          n !== r && (n = N(n)),
          (t = n),
          t
        );
      }
      function rt() {
        var t, n;
        return (
          (t = H),
          P.test(e.charAt(H))
            ? ((n = e.charAt(H)), H++)
            : ((n = r), W === 0 && q(ne)),
          n !== r && (n = f(n)),
          (t = n),
          t
        );
      }
      function it() {
        var t, n;
        return (
          (t = H),
          e.charCodeAt(H) === 42 ? ((n = F), H++) : ((n = r), W === 0 && q(I)),
          n !== r && (n = re(n)),
          (t = n),
          t
        );
      }
      function at() {
        var t, n, i, a;
        return (
          (t = H),
          (n = H),
          ie.test(e.charAt(H))
            ? ((i = e.charAt(H)), H++)
            : ((i = r), W === 0 && q(ae)),
          i === r
            ? ((H = n), (n = r))
            : (oe.test(e.charAt(H))
                ? ((a = e.charAt(H)), H++)
                : ((a = r), W === 0 && q(se)),
              a === r && (a = null),
              a === r ? ((H = n), (n = r)) : ((i = [i, a]), (n = i))),
          n !== r && (n = ce(n)),
          (t = n),
          t
        );
      }
      function ot() {
        var t, n, i, a, o;
        return (
          (t = H),
          (n = H),
          e.charCodeAt(H) === 37
            ? ((i = le), H++)
            : ((i = r), W === 0 && q(ue)),
          i === r
            ? ((H = n), (n = r))
            : (L.test(e.charAt(H))
                ? ((a = e.charAt(H)), H++)
                : ((a = r), W === 0 && q(R)),
              a === r
                ? ((H = n), (n = r))
                : (z.test(e.charAt(H))
                    ? ((o = e.charAt(H)), H++)
                    : ((o = r), W === 0 && q(B)),
                  o === r ? ((H = n), (n = r)) : ((i = [i, a, o]), (n = i)))),
          n === r &&
            (z.test(e.charAt(H))
              ? ((n = e.charAt(H)), H++)
              : ((n = r), W === 0 && q(B))),
          n !== r && (n = de(n)),
          (t = n),
          t
        );
      }
      function st() {
        var t, n, i, a, o, s, c;
        return (
          (t = H),
          (n = H),
          e.charCodeAt(H) === 64
            ? ((i = fe), H++)
            : ((i = r), W === 0 && q(pe)),
          i === r
            ? ((H = n), (n = r))
            : (e.charCodeAt(H) === 64
                ? ((a = fe), H++)
                : ((a = r), W === 0 && q(pe)),
              a === r &&
                ((a = H),
                e.substr(H, 2) === me
                  ? ((o = me), (H += 2))
                  : ((o = r), W === 0 && q(he)),
                o === r
                  ? ((H = a), (a = r))
                  : (ge.test(e.charAt(H))
                      ? ((s = e.charAt(H)), H++)
                      : ((s = r), W === 0 && q(V)),
                    s === r ? ((H = a), (a = r)) : ((o = [o, s]), (a = o))),
                a === r &&
                  ((a = H),
                  e.substr(H, 2) === _e
                    ? ((o = _e), (H += 2))
                    : ((o = r), W === 0 && q(ve)),
                  o === r
                    ? ((H = a), (a = r))
                    : (ge.test(e.charAt(H))
                        ? ((s = e.charAt(H)), H++)
                        : ((s = r), W === 0 && q(V)),
                      s === r ? ((H = a), (a = r)) : ((o = [o, s]), (a = o))),
                  a === r &&
                    ((a = H),
                    e.substr(H, 2) === ye
                      ? ((o = ye), (H += 2))
                      : ((o = r), W === 0 && q(be)),
                    o === r
                      ? ((H = a), (a = r))
                      : (xe.test(e.charAt(H))
                          ? ((s = e.charAt(H)), H++)
                          : ((s = r), W === 0 && q(Se)),
                        s === r ? ((H = a), (a = r)) : ((o = [o, s]), (a = o))),
                    a === r &&
                      ((a = H),
                      e.substr(H, 2) === Ce
                        ? ((o = Ce), (H += 2))
                        : ((o = r), W === 0 && q(we)),
                      o === r
                        ? ((H = a), (a = r))
                        : (L.test(e.charAt(H))
                            ? ((s = e.charAt(H)), H++)
                            : ((s = r), W === 0 && q(R)),
                          s === r
                            ? ((H = a), (a = r))
                            : (z.test(e.charAt(H))
                                ? ((c = e.charAt(H)), H++)
                                : ((c = r), W === 0 && q(B)),
                              c === r && (c = null),
                              c === r
                                ? ((H = a), (a = r))
                                : ((o = [o, s, c]), (a = o)))),
                      a === r &&
                        ((a = H),
                        e.substr(H, 2) === Te
                          ? ((o = Te), (H += 2))
                          : ((o = r), W === 0 && q(Ee)),
                        o === r
                          ? ((H = a), (a = r))
                          : (L.test(e.charAt(H))
                              ? ((s = e.charAt(H)), H++)
                              : ((s = r), W === 0 && q(R)),
                            s === r
                              ? ((H = a), (a = r))
                              : (z.test(e.charAt(H))
                                  ? ((c = e.charAt(H)), H++)
                                  : ((c = r), W === 0 && q(B)),
                                c === r && (c = null),
                                c === r
                                  ? ((H = a), (a = r))
                                  : ((o = [o, s, c]), (a = o))))))))),
              a === r && (a = null),
              a === r ? ((H = n), (n = r)) : ((i = [i, a]), (n = i))),
          n !== r && (n = De(n)),
          (t = n),
          t
        );
      }
      function ct() {
        var e, t;
        return (
          (e = H),
          (t = lt()),
          t === r && (t = ut()),
          t !== r && (t = Oe(t)),
          (e = t),
          e
        );
      }
      function lt() {
        var t, n, i, a, o, s;
        return (
          (t = H),
          (n = H),
          e.charCodeAt(H) === 43
            ? ((i = ke), H++)
            : ((i = r), W === 0 && q(Ae)),
          i === r
            ? ((H = n), (n = r))
            : (e.charCodeAt(H) === 43
                ? ((a = ke), H++)
                : ((a = r), W === 0 && q(Ae)),
              a === r &&
                ((a = H),
                L.test(e.charAt(H))
                  ? ((o = e.charAt(H)), H++)
                  : ((o = r), W === 0 && q(R)),
                o === r
                  ? ((H = a), (a = r))
                  : (z.test(e.charAt(H))
                      ? ((s = e.charAt(H)), H++)
                      : ((s = r), W === 0 && q(B)),
                    s === r && (s = null),
                    s === r ? ((H = a), (a = r)) : ((o = [o, s]), (a = o)))),
              a === r && (a = null),
              a === r ? ((H = n), (n = r)) : ((i = [i, a]), (n = i))),
          n !== r && (n = je(n)),
          (t = n),
          t
        );
      }
      function ut() {
        var t, n, i, a, o, s;
        return (
          (t = H),
          (n = H),
          e.charCodeAt(H) === 45
            ? ((i = Me), H++)
            : ((i = r), W === 0 && q(Ne)),
          i === r
            ? ((H = n), (n = r))
            : (e.charCodeAt(H) === 45
                ? ((a = Me), H++)
                : ((a = r), W === 0 && q(Ne)),
              a === r &&
                ((a = H),
                L.test(e.charAt(H))
                  ? ((o = e.charAt(H)), H++)
                  : ((o = r), W === 0 && q(R)),
                o === r
                  ? ((H = a), (a = r))
                  : (z.test(e.charAt(H))
                      ? ((s = e.charAt(H)), H++)
                      : ((s = r), W === 0 && q(B)),
                    s === r && (s = null),
                    s === r ? ((H = a), (a = r)) : ((o = [o, s]), (a = o)))),
              a === r && (a = null),
              a === r ? ((H = n), (n = r)) : ((i = [i, a]), (n = i))),
          n !== r && (n = Pe(n)),
          (t = n),
          t
        );
      }
      function dt() {
        var t, n, i, a;
        return (
          (t = H),
          (n = H),
          e.charCodeAt(H) === 72
            ? ((i = Fe), H++)
            : ((i = r), W === 0 && q(Ie)),
          i === r
            ? ((H = n), (n = r))
            : (z.test(e.charAt(H))
                ? ((a = e.charAt(H)), H++)
                : ((a = r), W === 0 && q(B)),
              a === r && (a = null),
              a === r ? ((H = n), (n = r)) : ((i = [i, a]), (n = i))),
          n !== r && (n = Le(n)),
          (t = n),
          t
        );
      }
      function ft() {
        var t, n, i, a, o, s, c;
        if (
          ((t = H),
          (n = H),
          e.charCodeAt(H) === 58
            ? ((i = Re), H++)
            : ((i = r), W === 0 && q(ze)),
          i !== r)
        ) {
          if (
            ((a = H),
            L.test(e.charAt(H))
              ? ((o = e.charAt(H)), H++)
              : ((o = r), W === 0 && q(R)),
            o !== r)
          ) {
            for (
              s = [],
                z.test(e.charAt(H))
                  ? ((c = e.charAt(H)), H++)
                  : ((c = r), W === 0 && q(B));
              c !== r;
            )
              (s.push(c),
                z.test(e.charAt(H))
                  ? ((c = e.charAt(H)), H++)
                  : ((c = r), W === 0 && q(B)));
            s === r ? ((H = a), (a = r)) : ((o = [o, s]), (a = o));
          } else ((H = a), (a = r));
          (a === r &&
            (Be.test(e.charAt(H))
              ? ((a = e.charAt(H)), H++)
              : ((a = r), W === 0 && q(Ve))),
            a === r ? ((H = n), (n = r)) : ((i = [i, a]), (n = i)));
        } else ((H = n), (n = r));
        return (n !== r && (n = He(n)), (t = n), t);
      }
      function pt() {
        var t, n, i, a, o;
        return (
          (t = H),
          (n = H),
          L.test(e.charAt(H))
            ? ((i = e.charAt(H)), H++)
            : ((i = r), W === 0 && q(R)),
          i === r
            ? ((H = n), (n = r))
            : (z.test(e.charAt(H))
                ? ((a = e.charAt(H)), H++)
                : ((a = r), W === 0 && q(B)),
              a === r && (a = null),
              a === r
                ? ((H = n), (n = r))
                : (z.test(e.charAt(H))
                    ? ((o = e.charAt(H)), H++)
                    : ((o = r), W === 0 && q(B)),
                  o === r && (o = null),
                  o === r ? ((H = n), (n = r)) : ((i = [i, a, o]), (n = i)))),
          n !== r && (n = Ue(n)),
          (t = n),
          t
        );
      }
      if (((Ke = a()), Ke !== r && H === e.length)) return Ke;
      throw (
        Ke !== r && H < e.length && q(qe()),
        Ze(
          Ge,
          U < e.length ? e.charAt(U) : null,
          U < e.length ? Ye(U, U + 1) : Ye(U, U),
        )
      );
    }
    return { SyntaxError: t, parse: n };
  })(),
  de = class {
    static parse(e) {
      return B.parse(e);
    }
  };
de.SyntaxError = B.SyntaxError;
var fe = {
    C2H4O2: `acetic acid`,
    C3H6O: `acetone`,
    C2H3N: `acetonitrile`,
    C6H6: `benzene`,
    CCl4: `carbon tetrachloride`,
    C6H5Cl: `chlorobenzene`,
    CHCl3: `chloroform`,
    C6H12: `cyclohexane`,
    C2H4Cl2: `1,2-dichloroethane`,
    C4H10O3: `diethylene glycol`,
    C6H14O3: `diglyme`,
    C4H10O2: `DME`,
    C3H7NO: `DMF`,
    C2H6OS: `DMSO`,
    C2H6O: `ethanol`,
    C2H6O2: `ethylene glycol`,
    C3H8O3: `glycerin`,
    C7H16: `heptane`,
    C6H18N3OP: `HMPA`,
    C6H18N3P: `HMPT`,
    C6H14: `hexane`,
    CH4O: `methanol`,
    C5H12O: `MTBE`,
    CH2Cl2: `methylene chloride`,
    CH5H9NO: `NMP`,
    CH3NO2: `nitromethane`,
    C5H12: `pentane`,
    C5H5N: `pyridine`,
    C7H8: `toluene`,
    C6H15N: `triethyl amine`,
    H2O: `water`,
  },
  pe = class {
    constructor(e, t) {
      ((this.drawer = new R(t)),
        (this.molOpts = this.drawer.opts),
        (this.defaultOptions = {
          scale: this.molOpts.scale > 0 ? this.molOpts.scale : 1,
          fontSize: this.molOpts.fontSizeLarge * 0.8,
          fontFamily: `Arial, Helvetica, sans-serif`,
          spacing: 10,
          plus: { size: 9, thickness: 1 },
          arrow: {
            length: this.molOpts.bondLength * 4,
            headSize: 6,
            thickness: 1,
            margin: 3,
          },
          weights: { normalize: !1 },
        }),
        (this.opts = O.extend(!0, this.defaultOptions, e)));
    }
    draw(e, t, n = `light`, r = null, i = `{reagents}`, a = ``, o = !1) {
      if (
        ((this.themeManager = new x(this.molOpts.themes, n)),
        this.opts.weights.normalize)
      ) {
        let e = -(2 ** 53 - 1),
          t = 2 ** 53 - 1;
        if (`reactants` in r)
          for (let n = 0; n < r.reactants.length; n++)
            for (let i = 0; i < r.reactants[n].length; i++)
              (r.reactants[n][i] < t && (t = r.reactants[n][i]),
                r.reactants[n][i] > e && (e = r.reactants[n][i]));
        if (`reagents` in r)
          for (let n = 0; n < r.reagents.length; n++)
            for (let i = 0; i < r.reagents[n].length; i++)
              (r.reagents[n][i] < t && (t = r.reagents[n][i]),
                r.reagents[n][i] > e && (e = r.reagents[n][i]));
        if (`products` in r)
          for (let n = 0; n < r.products.length; n++)
            for (let i = 0; i < r.products[n].length; i++)
              (r.products[n][i] < t && (t = r.products[n][i]),
                r.products[n][i] > e && (e = r.products[n][i]));
        let n = Math.max(Math.abs(t), Math.abs(e));
        if ((n === 0 && (n = 1), `reactants` in r))
          for (let e = 0; e < r.reactants.length; e++)
            for (let t = 0; t < r.reactants[e].length; t++)
              r.reactants[e][t] /= n;
        if (`reagents` in r)
          for (let e = 0; e < r.reagents.length; e++)
            for (let t = 0; t < r.reagents[e].length; t++)
              r.reagents[e][t] /= n;
        if (`products` in r)
          for (let e = 0; e < r.products.length; e++)
            for (let t = 0; t < r.products[e].length; t++)
              r.products[e][t] /= n;
      }
      let s = null;
      for (
        t === null || t === `svg`
          ? ((s = document.createElementNS(
              `http://www.w3.org/2000/svg`,
              `svg`,
            )),
            s.setAttribute(`xmlns`, `http://www.w3.org/2000/svg`),
            s.setAttributeNS(null, `width`, `500`),
            s.setAttributeNS(null, `height`, `500`))
          : (s =
              typeof t == `string` || t instanceof String
                ? document.getElementById(t)
                : t);
        s.firstChild;
      )
        s.removeChild(s.firstChild);
      let c = [],
        l = 0;
      for (let t = 0; t < e.reactants.length; t++) {
        t > 0 &&
          c.push({
            width: this.opts.plus.size * this.opts.scale,
            height: this.opts.plus.size * this.opts.scale,
            svg: this.getPlus(),
          });
        let i = null;
        r && `reactants` in r && r.reactants.length > t && (i = r.reactants[t]);
        let a = document.createElementNS(`http://www.w3.org/2000/svg`, `svg`);
        this.drawer.draw(
          e.reactants[t],
          a,
          n,
          i,
          o,
          [],
          this.opts.weights.normalize,
        );
        let s = {
          width: a.viewBox.baseVal.width * this.opts.scale,
          height: a.viewBox.baseVal.height * this.opts.scale,
          svg: a,
        };
        (c.push(s), s.height > l && (l = s.height));
      }
      c.push({
        width: this.opts.arrow.length * this.opts.scale,
        height: this.opts.arrow.headSize * this.opts.scale * 2,
        svg: this.getArrow(),
      });
      let u = ``;
      for (let t = 0; t < e.reagents.length; t++) {
        t > 0 && (u += `, `);
        let n = this.drawer.getMolecularFormula(e.reagents[t]);
        (n in fe && (n = fe[n]), (u += L.replaceNumbersWithSubscript(n)));
      }
      i = i.replace(`{reagents}`, u);
      let d = L.writeText(
          i,
          this.themeManager,
          this.opts.fontSize * this.opts.scale,
          this.opts.fontFamily,
          this.opts.arrow.length * this.opts.scale,
        ),
        f = (this.opts.arrow.length * this.opts.scale - d.width) / 2;
      c.push({
        svg: d.svg,
        height: d.height,
        width: d.width,
        offsetX:
          -(this.opts.arrow.length * this.opts.scale + this.opts.spacing) + f,
        offsetY: -(d.height / 2) - this.opts.arrow.margin,
        position: `relative`,
      });
      let p = L.writeText(
        a,
        this.themeManager,
        this.opts.fontSize * this.opts.scale,
        this.opts.fontFamily,
        this.opts.arrow.length * this.opts.scale,
      );
      ((f = (this.opts.arrow.length * this.opts.scale - p.width) / 2),
        c.push({
          svg: p.svg,
          height: p.height,
          width: p.width,
          offsetX:
            -(this.opts.arrow.length * this.opts.scale + this.opts.spacing) + f,
          offsetY: p.height / 2 + this.opts.arrow.margin,
          position: `relative`,
        }));
      for (let t = 0; t < e.products.length; t++) {
        t > 0 &&
          c.push({
            width: this.opts.plus.size * this.opts.scale,
            height: this.opts.plus.size * this.opts.scale,
            svg: this.getPlus(),
          });
        let i = null;
        r && `products` in r && r.products.length > t && (i = r.products[t]);
        let a = document.createElementNS(`http://www.w3.org/2000/svg`, `svg`);
        this.drawer.draw(
          e.products[t],
          a,
          n,
          i,
          o,
          [],
          this.opts.weights.normalize,
        );
        let s = {
          width: a.viewBox.baseVal.width * this.opts.scale,
          height: a.viewBox.baseVal.height * this.opts.scale,
          svg: a,
        };
        (c.push(s), s.height > l && (l = s.height));
      }
      let m = 0,
        h = 0,
        g = 0;
      c.forEach((e) => {
        let t = e.offsetX || 0,
          n = e.offsetY || 0,
          r = (l - e.height) / 2 + n;
        ((h = Math.max(h, r + e.height)),
          (m = Math.min(m, r)),
          e.svg.setAttributeNS(null, `x`, Math.round(g + t)),
          e.svg.setAttributeNS(null, `y`, Math.round(r)),
          e.svg.setAttributeNS(null, `width`, Math.round(e.width)),
          e.svg.setAttributeNS(null, `height`, Math.round(e.height)),
          s.appendChild(e.svg),
          e.position !== `relative` &&
            (g += Math.round(e.width + this.opts.spacing + t)));
      });
      let _ = Math.max(l, h - m);
      return (
        s.setAttributeNS(null, `viewBox`, `0 ${m} ${g} ${_}`),
        (s.style.width = g + `px`),
        (s.style.height = l + `px`),
        s
      );
    }
    getPlus() {
      let e = this.opts.plus.size,
        t = this.opts.plus.thickness,
        n = document.createElementNS(`http://www.w3.org/2000/svg`, `svg`),
        r = document.createElementNS(`http://www.w3.org/2000/svg`, `rect`),
        i = document.createElementNS(`http://www.w3.org/2000/svg`, `rect`);
      return (
        n.setAttributeNS(null, `id`, `plus`),
        r.setAttributeNS(null, `x`, 0),
        r.setAttributeNS(null, `y`, e / 2 - t / 2),
        r.setAttributeNS(null, `width`, e),
        r.setAttributeNS(null, `height`, t),
        r.setAttributeNS(null, `fill`, this.themeManager.getColor(`C`)),
        i.setAttributeNS(null, `x`, e / 2 - t / 2),
        i.setAttributeNS(null, `y`, 0),
        i.setAttributeNS(null, `width`, t),
        i.setAttributeNS(null, `height`, e),
        i.setAttributeNS(null, `fill`, this.themeManager.getColor(`C`)),
        n.appendChild(r),
        n.appendChild(i),
        n.setAttributeNS(null, `viewBox`, `0 0 ${e} ${e}`),
        n
      );
    }
    getArrowhead() {
      let e = this.opts.arrow.headSize,
        t = document.createElementNS(`http://www.w3.org/2000/svg`, `marker`),
        n = document.createElementNS(`http://www.w3.org/2000/svg`, `polygon`);
      return (
        t.setAttributeNS(null, `id`, `arrowhead`),
        t.setAttributeNS(null, `viewBox`, `0 0 ${e} ${e}`),
        t.setAttributeNS(null, `markerUnits`, `userSpaceOnUse`),
        t.setAttributeNS(null, `markerWidth`, e),
        t.setAttributeNS(null, `markerHeight`, e),
        t.setAttributeNS(null, `refX`, 0),
        t.setAttributeNS(null, `refY`, e / 2),
        t.setAttributeNS(null, `orient`, `auto`),
        t.setAttributeNS(null, `fill`, this.themeManager.getColor(`C`)),
        n.setAttributeNS(null, `points`, `0 0, ${e} ${e / 2}, 0 ${e}`),
        t.appendChild(n),
        t
      );
    }
    getCDArrowhead() {
      let e = this.opts.arrow.headSize,
        t = (7 / 4.5) * e,
        n = document.createElementNS(`http://www.w3.org/2000/svg`, `marker`),
        r = document.createElementNS(`http://www.w3.org/2000/svg`, `path`);
      return (
        n.setAttributeNS(null, `id`, `arrowhead`),
        n.setAttributeNS(null, `viewBox`, `0 0 ${t} ${e}`),
        n.setAttributeNS(null, `markerUnits`, `userSpaceOnUse`),
        n.setAttributeNS(null, `markerWidth`, t * 2),
        n.setAttributeNS(null, `markerHeight`, e * 2),
        n.setAttributeNS(null, `refX`, 2.2),
        n.setAttributeNS(null, `refY`, 2.2),
        n.setAttributeNS(null, `orient`, `auto`),
        n.setAttributeNS(null, `fill`, this.themeManager.getColor(`C`)),
        r.setAttributeNS(null, `style`, `fill-rule:nonzero;`),
        r.setAttributeNS(
          null,
          `d`,
          `m 0 0 l 7 2.25 l -7 2.25 c 0 0 0.735 -1.084 0.735 -2.28 c 0 -1.196 -0.735 -2.22 -0.735 -2.22 z`,
        ),
        n.appendChild(r),
        n
      );
    }
    getArrow() {
      let e = this.opts.arrow.headSize,
        t = this.opts.arrow.length,
        n = document.createElementNS(`http://www.w3.org/2000/svg`, `svg`),
        r = document.createElementNS(`http://www.w3.org/2000/svg`, `defs`),
        i = document.createElementNS(`http://www.w3.org/2000/svg`, `line`);
      return (
        r.appendChild(this.getCDArrowhead()),
        n.appendChild(r),
        n.setAttributeNS(null, `id`, `arrow`),
        i.setAttributeNS(null, `x1`, 0),
        i.setAttributeNS(null, `y1`, -this.opts.arrow.thickness / 2),
        i.setAttributeNS(null, `x2`, t),
        i.setAttributeNS(null, `y2`, -this.opts.arrow.thickness / 2),
        i.setAttributeNS(null, `stroke-width`, this.opts.arrow.thickness),
        i.setAttributeNS(null, `stroke`, this.themeManager.getColor(`C`)),
        i.setAttributeNS(null, `marker-end`, `url(#arrowhead)`),
        n.appendChild(i),
        n.setAttributeNS(
          null,
          `viewBox`,
          `0 ${-e / 2} ${t + (7 / 4.5) * e} ${e}`,
        ),
        n
      );
    }
  },
  me = class {
    constructor(e) {
      ((this.reactantsSmiles = []),
        (this.reagentsSmiles = []),
        (this.productsSmiles = []),
        (this.reactantsWeights = []),
        (this.reagentsWeights = []),
        (this.productsWeights = []),
        (this.reactants = []),
        (this.reagents = []),
        (this.products = []));
      let t = e.split(`>`);
      if (t.length !== 3)
        throw Error(
          `Invalid reaction SMILES: Expected exactly two ">" symbols.`,
        );
      (t[0] !== `` && (this.reactantsSmiles = t[0].split(`.`)),
        t[1] !== `` && (this.reagentsSmiles = t[1].split(`.`)),
        t[2] !== `` && (this.productsSmiles = t[2].split(`.`)));
      for (let e = 0; e < this.reactantsSmiles.length; e++)
        this.reactants.push(B.parse(this.reactantsSmiles[e]));
      for (let e = 0; e < this.reagentsSmiles.length; e++)
        this.reagents.push(B.parse(this.reagentsSmiles[e]));
      for (let e = 0; e < this.productsSmiles.length; e++)
        this.products.push(B.parse(this.productsSmiles[e]));
    }
  },
  he = class {
    static parse(e) {
      return new me(e);
    }
  },
  ge = class e {
    constructor(e = {}, t = {}) {
      ((this.drawer = new R(e)),
        (this.reactionDrawer = new pe(
          t,
          JSON.parse(JSON.stringify(this.drawer.opts)),
        )));
    }
    static apply(
      t = {},
      n = {},
      r = `data-smiles`,
      i = `light`,
      a = null,
      o = null,
    ) {
      new e(t, n).apply(r, i, a, o);
    }
    apply(t = `data-smiles`, n = `light`, r = null, i = null) {
      document.querySelectorAll(`[${t}]`).forEach((a) => {
        let o = a.getAttribute(t);
        if (o === null) throw Error(`No SMILES provided.`);
        let s = n,
          c = null;
        if (
          (a.hasAttribute(`data-smiles-theme`) &&
            (s = a.getAttribute(`data-smiles-theme`)),
          a.hasAttribute(`data-smiles-weights`) &&
            (c = a
              .getAttribute(`data-smiles-weights`)
              .split(`,`)
              .map(parseFloat)),
          (a.hasAttribute(`data-smiles-reactant-weights`) ||
            a.hasAttribute(`data-smiles-reagent-weights`) ||
            a.hasAttribute(`data-smiles-product-weights`)) &&
            ((c = { reactants: [], reagents: [], products: [] }),
            a.hasAttribute(`data-smiles-reactant-weights`) &&
              (c.reactants = a
                .getAttribute(`data-smiles-reactant-weights`)
                .split(`;`)
                .map((e) => e.split(`,`).map(parseFloat))),
            a.hasAttribute(`data-smiles-reagent-weights`) &&
              (c.reagents = a
                .getAttribute(`data-smiles-reagent-weights`)
                .split(`;`)
                .map((e) => e.split(`,`).map(parseFloat))),
            a.hasAttribute(`data-smiles-product-weights`) &&
              (c.products = a
                .getAttribute(`data-smiles-product-weights`)
                .split(`;`)
                .map((e) => e.split(`,`).map(parseFloat)))),
          a.hasAttribute(`data-smiles-options`) ||
            a.hasAttribute(`data-smiles-reaction-options`))
        ) {
          let t = {};
          if (a.hasAttribute(`data-smiles-options`)) {
            let e = a.getAttribute(`data-smiles-options`);
            try {
              t = JSON.parse(e);
            } catch {
              t = JSON.parse(e.replace(/'/g, `"`));
            }
          }
          let n = {};
          if (a.hasAttribute(`data-smiles-reaction-options`)) {
            let e = a.getAttribute(`data-smiles-reaction-options`);
            try {
              n = JSON.parse(e);
            } catch {
              n = JSON.parse(e.replace(/'/g, `"`));
            }
          }
          new e(t, n).draw(o, a, s, r, i, c);
        } else this.draw(o, a, s, r, i, c);
      });
    }
    draw(e, t, n = `light`, r = null, i = null, a = null) {
      let o = [];
      [e, ...o] = e.split(` `);
      let s = o.join(` `),
        c = {};
      if (s.includes(`__`)) {
        let e = s.substring(s.indexOf(`__`) + 2, s.lastIndexOf(`__`));
        try {
          c = JSON.parse(e);
        } catch {
          c = JSON.parse(e.replace(/'/g, `"`));
        }
      }
      if (
        ((c = O.extend(
          !0,
          { textAboveArrow: `{reagents}`, textBelowArrow: `` },
          c,
        )),
        e.includes(`>`))
      )
        try {
          this.drawReaction(e, t, n, c, a, r);
        } catch (e) {
          i ? i(e) : console.error(e);
        }
      else
        try {
          this.drawMolecule(e, t, n, a, r);
        } catch (e) {
          i ? i(e) : console.error(e);
        }
    }
    drawMolecule(e, t, n, r, i) {
      let a = B.parse(e);
      if (t === null || t === `svg`) {
        let e = this.drawer.draw(a, null, n, r),
          t = this.getDimensions(e);
        (e.setAttributeNS(null, `width`, `` + t.w),
          e.setAttributeNS(null, `height`, `` + t.h),
          i && i(e));
      } else if (t === `canvas`) {
        let e = this.svgToCanvas(this.drawer.draw(a, null, n, r));
        i && i(e);
      } else if (t === `img`) {
        let e = this.svgToImg(this.drawer.draw(a, null, n, r));
        i && i(e);
      } else
        t instanceof HTMLImageElement
          ? (this.svgToImg(this.drawer.draw(a, null, n, r), t), i && i(t))
          : t instanceof SVGElement
            ? (this.drawer.draw(a, t, n, r), i && i(t))
            : document.querySelectorAll(t).forEach((e) => {
                let t = e.nodeName.toLowerCase();
                t === `svg`
                  ? (this.drawer.draw(a, e, n, r), i && i(e))
                  : t === `canvas`
                    ? (this.svgToCanvas(this.drawer.draw(a, null, n, r), e),
                      i && i(e))
                    : t === `img` &&
                      (this.svgToImg(this.drawer.draw(a, null, n, r), e),
                      i && i(e));
              });
    }
    drawReaction(e, t, n, r, i, a) {
      let o = he.parse(e);
      if (t === null || t === `svg`) {
        let e = this.reactionDrawer.draw(o, null, n),
          t = this.getDimensions(e);
        (e.setAttributeNS(null, `width`, `` + t.w),
          e.setAttributeNS(null, `height`, `` + t.h),
          a && a(e));
      } else if (t === `canvas`) {
        let e = this.svgToCanvas(
          this.reactionDrawer.draw(
            o,
            null,
            n,
            i,
            r.textAboveArrow,
            r.textBelowArrow,
          ),
        );
        a && a(e);
      } else if (t === `img`) {
        let e = this.svgToImg(
          this.reactionDrawer.draw(
            o,
            null,
            n,
            i,
            r.textAboveArrow,
            r.textBelowArrow,
          ),
        );
        a && a(e);
      } else
        t instanceof HTMLImageElement
          ? (this.svgToImg(
              this.reactionDrawer.draw(
                o,
                null,
                n,
                i,
                r.textAboveArrow,
                r.textBelowArrow,
              ),
              t,
            ),
            a && a(t))
          : t instanceof SVGElement
            ? (this.reactionDrawer.draw(
                o,
                t,
                n,
                i,
                r.textAboveArrow,
                r.textBelowArrow,
              ),
              a && a(t))
            : document.querySelectorAll(t).forEach((e) => {
                let t = e.nodeName.toLowerCase();
                t === `svg`
                  ? (this.reactionDrawer.draw(
                      o,
                      e,
                      n,
                      i,
                      r.textAboveArrow,
                      r.textBelowArrow,
                    ),
                    this.reactionDrawer.opts.scale <= 0 &&
                      ((e.style.width = null), (e.style.height = null)),
                    a && a(e))
                  : t === `canvas`
                    ? (this.svgToCanvas(
                        this.reactionDrawer.draw(
                          o,
                          null,
                          n,
                          i,
                          r.textAboveArrow,
                          r.textBelowArrow,
                        ),
                        e,
                      ),
                      a && a(e))
                    : t === `img` &&
                      (this.svgToImg(
                        this.reactionDrawer.draw(
                          o,
                          null,
                          n,
                          i,
                          r.textAboveArrow,
                          r.textBelowArrow,
                        ),
                        e,
                      ),
                      a && a(e));
              });
    }
    svgToCanvas(e, t = null) {
      t === null && (t = document.createElement(`canvas`));
      let n = this.getDimensions(t, e);
      return (L.svgToCanvas(e, t, n.w, n.h), t);
    }
    svgToImg(e, t = null) {
      t === null && (t = document.createElement(`img`));
      let n = this.getDimensions(t, e);
      return (L.svgToImg(e, t, n.w, n.h), t);
    }
    getDimensions(e, t = null) {
      let n = this.drawer.opts.width,
        r = this.drawer.opts.height;
      return (
        this.drawer.opts.scale <= 0
          ? (e instanceof SVGElement ||
              (n === null && (n = e.width), r === null && (r = e.height)),
            e.style.width !== `` && (n = parseInt(e.style.width)),
            e.style.height !== `` && (r = parseInt(e.style.height)))
          : t &&
            ((n = parseFloat(t.style.width)), (r = parseFloat(t.style.height))),
        { w: n, h: r }
      );
    }
  },
  V = class e {
    static clean(e) {
      return e.replace(/[^A-Za-z0-9@.+\-?!()[\]{}/\\=#$:*]/g, ``);
    }
    static apply(t, n = `canvas[data-smiles]`, r = `light`, i) {
      let a = new z(t),
        o = document.querySelectorAll(n);
      for (let t = 0; t < o.length; t++) {
        let n = o[t];
        e.parse(
          n.getAttribute(`data-smiles`),
          function (e) {
            a.draw(e, n, r, !1);
          },
          function (e) {
            i && i(e);
          },
        );
      }
    }
    static parse(e, t, n) {
      try {
        t && t(de.parse(e));
      } catch (e) {
        n && n(e);
      }
    }
    static parseReaction(e, t, n) {
      try {
        t && t(he.parse(e));
      } catch (e) {
        n && n(e);
      }
    }
  };
((V.Version = `2.4.1`),
  (V.Drawer = z),
  (V.GaussDrawer = le),
  (V.Parser = de),
  (V.ReactionDrawer = pe),
  (V.ReactionParser = he),
  (V.SmiDrawer = ge),
  (V.SvgDrawer = R));
var _e = V;
typeof window < `u` &&
  window.document &&
  window.document.createElement &&
  ((window.SmilesDrawer = _e), (window.SmiDrawer = ge));
export { _e as default };
