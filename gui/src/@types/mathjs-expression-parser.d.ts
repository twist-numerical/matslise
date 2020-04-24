declare let MEP: {
  eval(text: string, scope?: any): number;
  compile(
    text: string
  ): {
    eval: (scope?: any) => number;
  };
};

declare module "mathjs-expression-parser" {
  export = MEP;
}
