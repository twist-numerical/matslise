import MathParser from './MathParser';

Math.square = (x) => x*x;

const DF = function(f, df, toString) {
	const r = function(v) {
		return DF.chain(r, v);
	};
	r.get = f;
	r.diff = df;
	r.toString = toString;
	r.toFunction = function (...args) {
		// eslint-disable-next-line 
		return new Function(...args, "return "+this.toString());
	};
	r.toValue = function() {
		return r.toFunction()();
	}
	return r;
}

DF.number = (n) => DF(() => n, () => DF.zero, () => ""+n);
const N0 = DF.zero = DF.number(0);
const N1 = DF.one = DF.number(1);
const X = DF.x = DF((subs, x) => x, () => N1, (s) => s)
DF.variable = (name) => DF((subs, x) => subs[name], (v) => v === name?1:0, (s) => name);
DF.chain = (f, g) => DF(
	(subs, x) => f.get(subs, g.get(subs, x)),
	v => DF.multiply(DF.chain(f.diff(v), g), g.diff(v)),
	(s) => f.toString(g.toString(s)));
DF.negate = (f) => f === N0?N0:DF((subs, x) => -f.get(subs, x), (v) => DF.negate(f.diff(v)), (s) => "(-"+f.toString(s)+")");
DF.sin = DF((subs, x) => Math.sin(x), () => DF.cos, (s) => "Math.sin("+s+")");
DF.cos = DF((subs, x) => Math.cos(x), () => DF.negate(DF.sin), (s) => "Math.cos("+s+")");
DF.tan = DF((subs, x) => Math.tan(x), () => DF.divide(N1, DF.square(DF.cos)), (s) => "Math.tan("+s+")");
DF.atan = DF((subs, x) => Math.atan(x), () => DF.divide(N1, DF.add(N1, DF.square(X))), (s) => "Math.atan("+s+")");
DF.exp = DF((subs, x) => Math.exp(x), () => DF.exp, (s) => "Math.exp("+s+")");
DF.multiply = (f, g) => f === N0 || g === N0?N0:f===N1?g:g===N1?f:f===g?DF.square(f):DF(
	(subs, x) => f.get(subs, x) * g.get(subs, x),
	v => DF.add(DF.multiply(f.diff(v),g), DF.multiply(f, g.diff(v))),
	(s) => "("+f.toString(s)+"*"+g.toString(s)+")");
DF.divide = (f, g) => f === N0?N0:g===N1?f:DF(
	(subs, x) => f.get(subs, x) / g.get(subs, x),
	v => DF.divide(DF.subtract(DF.multiply(f.diff(v),g), DF.multiply(f, g.diff(v))), DF.square(g)),
	(s) => "("+f.toString(s)+"/"+g.toString(s)+")");
DF.subtract = (f, g) => f===N0?DF.negate(g):g===N0?f:f === g?N0:DF(
	(subs, x) => f.get(subs, x) - g.get(subs, x),
	v => DF.subtract(f.diff(v), g.diff(v)),
	(s) => "("+f.toString(s)+" - "+g.toString(s)+")");
DF.add = (f, g) => f === N0?g:g === N0?f:DF(
	(subs, x) => f.get(subs, x).add(g.get(subs, x)),
	v => DF.add(f.diff(v), g.diff(v)),
	(s) => "("+f.toString(s)+" + "+g.toString(s)+")");
DF.square = (f) => DF(
	(subs, x) => Math.square(f.get(subs, x)),
	v => DF.multiply(DF.number(2), DF.multiply(f.diff(v), f)),
	(s) => "Math.square("+f.toString(s)+")");
DF.pow = (f, e) => DF(
	(subs, x) => Math.pow(f.get(subs, x), e),
	v => DF.multiply(e, DF.multiply(f.diff(v), DF.pow(f, e-1))),
	(s) => "Math.pow("+f.toString(s)+","+e+")");
DF.sqrt = (f) => {
	let pow = DF.pow(f, .5);
	pow.toString = (s) => "Math.sqrt("+f.toString(s)+")";
	return pow; 
};
DF.ifelse = (crit, f, g) => DF(
	(subs, x) => crit.get(subs, x) >= 0? f.get(subs,x) : g.get(subs,x),
	v => DF.ifelse(crit, f.diff(v), g.diff(v)),
	(s) => "(("+crit.toString(s)+" >= 0) ? "+f.toString(s)+" : "+g.toString(s)+")"
	);
DF.identity = (f) => f;

class DFParser extends MathParser {
	constants = {
		pi: DF.number(Math.PI),
	}
	constructor() {
		super(DF.number,
			v => this.constants[v],
			(() => {
				const operators = {};
				["add", "subtract", "multiply", "divide", "negate", "pow", "identity"].forEach(
					f => operators[f] = DF[f]);
				return operators;
			})(),
			(() => {
				const functions = {};
				["sqrt", "sin", "cos", "tan", "atan", "exp", "square"].forEach(
					f => functions[f] = DF[f]);
				return functions;
			})());
	}

	parse(value, variables = []) {
		this.variables = (v) => {
			if(variables.indexOf(v) >= 0)
				return DF.variable(v);
			return this.constants[v];
		};
		const parsed = super.parse(value);
		parsed.value = value;
		return parsed;
	}
}

const parser = new DFParser();
DF.parse = parser.parse.bind(parser);

export default DF;