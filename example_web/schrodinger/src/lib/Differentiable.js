import getParser from './getParser';

Number.prototype.get = function () { 
	return this; 
};
Number.prototype.diff = () => 0;
Number.prototype.toFunction = () => () => 0;
Array.prototype.get = function(subs, x) {
	return this.map((v) => v.get(subs, x)); 
};
Array.prototype.diff = function(v) { 
	return this.map((x) => x.diff(v));
};
Array.prototype.toFunction = function(...args) { 
	return new Function(...args, "return "+this.toString());
};

const DF = function(f, df, toString) {
	const r = function(v) {
		return DF.chain(r, v);
	};
	r.get = f;
	r.diff = df;
	r.toString = toString;
	r.toFunction = function (...args) {
		return new Function(...args, "return "+this.toString());
	};
	return r;
}

DF.variable = (name) => name === "PI" ? Math.PI : DF((subs, x) => subs[name], (v) => v === name?1:0, (s) => name);
DF.chain = (f, g) => DF(
	(subs, x) => f.get(subs, g.get(subs, x)),
	v => DF.multiply(DF.chain(f.diff(v), g), g.diff(v)),
	(s) => f.toString(g.toString(s)));

DF.negate = (f) => f === 0?0:DF((subs, x) => -f.get(subs, x), (v) => DF.negate(f.diff(v)), (s) => "(-"+f.toString(s)+")");
DF.sin = DF((subs, x) => Math.sin(x), () => DF.cos, (s) => "Math.sin("+s+")");
DF.cos = DF((subs, x) => Math.cos(x), () => DF.negate(DF.sin), (s) => "Math.cos("+s+")");
DF.exp = DF((subs, x) => Math.exp(x), () => DF.exp, (s) => "Math.exp("+s+")");
DF.multiply = (f, g) => f === 0 || g === 0?0:f===1?g:g===1?f:f===g?DF.square(f):DF(
	(subs, x) => f.get(subs, x) * g.get(subs, x),
	v => DF.add(DF.multiply(f.diff(v),g), DF.multiply(f, g.diff(v))),
	(s) => "("+f.toString(s)+"*"+g.toString(s)+")");
DF.divide = (f, g) => f === 0?0:g===1?f:DF(
	(subs, x) => f.get(subs, x) / g.get(subs, x),
	v => DF.divide(DF.subtract(DF.multiply(f.diff(v),g), DF.multiply(f, g.diff(v))), DF.square(g)),
	(s) => "("+f.toString(s)+"/"+g.toString(s)+")");
DF.subtract = (f, g) => f===0?DF.negate(g):g===0?f:f === g?0:DF(
	(subs, x) => f.get(subs, x) - g.get(subs, x),
	v => DF.subtract(f.diff(v), g.diff(v)),
	(s) => "("+f.toString(s)+" - "+g.toString(s)+")");
DF.add = (f, g) => f === 0?g:g === 0?f:DF(
	(subs, x) => f.get(subs, x).add(g.get(subs, x)),
	v => DF.add(f.diff(v), g.diff(v)),
	(s) => "("+f.toString(s)+" + "+g.toString(s)+")");
DF.square = (f) => DF(
	(subs, x) => Math.square(f.get(subs, x)),
	v => DF.multiply(2, DF.multiply(f.diff(v), f)),
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

let parser = getParser(
	(() => {
		const vars = {};
		return (v) => {
			if(!vars[v])
				vars[v] = DF.variable(v);
			return vars[v];
		};
	})(),
	(() => {
		const operators = {};
		["add", "subtract", "multiply", "divide", "negate", "pow", "identity"].forEach(
			f => operators[f] = DF[f]);
		return operators;
	})(),
	(() => {
		const functions = {};
		["sqrt", "sin", "cos", "exp", "square"].forEach(
			f => functions[f] = DF[f]);
		return functions;
	})());
DF.parse = parser.parse.bind(parser);

export default DF;