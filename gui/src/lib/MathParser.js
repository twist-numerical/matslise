class Parser {
	constructor(f) {
		this.f = f;
	}

	apply(s) {
		return this.f(s);
	}

	parse(str) {
		const applied = Parser.whitespaceWrap(this).apply(str);
		const arr = applied.filter(([v, s]) => s.length === 0).map(([v,s]) => v);
		if(arr.length === 1)
			return arr[0];

		if(arr.length === 0) {
			const error = new Error("Could not parse the given query");
			error.todo = str;
			applied.forEach(([_, s]) => {
				if(s.length < error.todo.length)
					error.todo = s;
			});
			error.done = str.substr(0, str.length - error.todo.length);
			throw error;
		}
		else
			throw new Error("The query is not uniquely defined");
	}

	bind(ap) {
		return new Parser(s => {
			let arr = [];
			let app = this.apply(s);
			app.forEach(([v, s]) => {
				arr = arr.concat(ap(v).apply(s));
			});
			return arr;
		})
	}

	chain(p) {
		return new Parser(s => {
			let arr = [];
			this.apply(s).forEach(([v, s]) => {
				arr = arr.concat(p.apply(s));
			});
			return arr;
		})
	}

	mplus(fp) {
		return new Parser(s => {
			return this.apply(s).concat(fp.apply(s));
		});
	}
}

let toString = (f) => { return {toString: () => f}; };

Parser.mplus = (a, ...all) => all.reduce((c,b) => c.mplus(b), a);

Parser.empty = new Parser(s => []);
Parser.pure = v => new Parser(s => [[v, s]]);

Parser.match = m => new Parser(s => {
	if(s.substr(0, m.length) === m)
		return [[m, s.substr(m.length)]];
	else
		return [];
});

Parser.plus = p => p.bind(v =>
	Parser.pure([v]).mplus(
		Parser.plus(p).bind(w =>
			Parser.pure([v].concat(w)))));
Parser.star = p => Parser.pure([]).mplus(Parser.plus(p));

Parser.regex = regex => new Parser(s => {
	let found = s.match(regex);
	if(found == null)
		return [];
	else
		return [[found[0], s.substr(found[0].length)]];
});

Parser.unsignedFloat
= Parser.regex(/^(([0-9]+(\.[0-9]*)?)|([0-9]*\.[0-9]+))(e-?[0-9]+)?/g);

Parser.word
= Parser.regex(/^[a-zA-Z][a-zA-Z0-9_]*/g);

Parser.whitespace = Parser.regex(/^\s*/g);

Parser.whitespaceWrap = (p) => Parser.whitespace.chain(
	p.bind(v =>
		Parser.whitespace.chain(
			Parser.pure(v))));

class MathParser {
	numbers = null;
	variables = null;
	operators = null;
	functions = null;

	constructor(
		numbers = (n) => ""+n,
		variables = (name) => toString(name),
		operators = {
			identity: (a) => a,
			negate: (a) => toString("-"+a),
			add: (a,b) => toString("("+a+" + "+b+")"),
			multiply: (a,b) => toString("("+a+" * "+b+")"),
			subtract: (a,b) => toString("("+a+" - "+b+")"),
			divide: (a,b) => toString("("+a+" / "+b+")"),
			pow: (a,b) => toString("("+a+"^"+b+")"),
		},
		functions = {}) {
		operators.powNegate = (a, b) => operators.pow(a, operators.negate(b));
		operators.powIdentity = (a, b) => operators.pow(a, operators.identity(b));

		this.numbers = numbers;
		this.variables = variables;
		this.operators = operators;
		this.functions = functions;
		this.parser = this.buildParser();
	}

	buildParser() {
		const cumul = [];

		const unaryExp = (priority, ops) => ops.reduce((c, [op, f]) =>
			c.mplus(Parser.regex(op).bind(_ =>
				cumul[priority-1].bind(v =>
					Parser.pure(this.operators[f](v))))),
			Parser.empty);

		const binExp = (priority, ops) => {
			let opParser = ops.reduce((c, [op, f]) =>
				c.mplus(Parser.whitespaceWrap(Parser.regex(op)).chain(Parser.pure(f))),
				Parser.empty);
			return Parser.pure(null).bind(_ => cumul[priority-1].bind(first =>
				Parser.plus(opParser.bind(f => cumul[priority-1]
					.bind(d => Parser.pure([f, d])))
				).bind(st => Parser.pure(
					st.reduce((c, [f, d]) => this.operators[f](c, d), first)))));
		};

		const parsers = [
		Parser.mplus(
			Parser.unsignedFloat.bind(n => Parser.pure(this.numbers(n))),
			Parser.word.bind(f => Parser.regex(/^\s*\(\s*/g).chain(
				cumul[cumul.length-1].bind(e => Parser.regex(/^\s*\)/g).bind(_=> {
					if (this.functions[f] === undefined)
						return Parser.empty;
					return Parser.pure(this.functions[f](e))
				})))),
			Parser.word.bind(v => {
				const parsedVar = this.variables(v);
				if(parsedVar === undefined)
					return Parser.empty;
				return Parser.pure(parsedVar);
			}),
			Parser.regex(/^\(\s*/g).bind(_ =>
				cumul[cumul.length-1].bind(v =>
					Parser.regex(/^\s*\)/g).chain(
						Parser.pure(v))))),
		null,
		binExp(2, [[/^\^/, "pow"]]),
		binExp(3, [[/^\^\s*-/, "powNegate"], [/^\^\s*\+/, "powIdentity"]]),
		unaryExp(4, [[/^-/, "negate"], [/^\+/, "identity"]]),
		binExp(5, [[/^\*/, "multiply"], [/^\//, "divide"]]),
		binExp(6, [[/^\+/, "add"], [/^-/, "subtract"]])
		];


		cumul.push(parsers[0]);
		for(let i = 1; i < parsers.length; ++i) {
			if(parsers[i] == null)
				cumul.push(cumul[i-1]);
			else
				cumul.push(cumul[i-1].mplus(parsers[i]));
		}

		return cumul[cumul.length-1];
	};

	parse(value) {
		return this.parser.parse(value);
	}
};
MathParser.Parser = Parser;

export default MathParser;
