import { n as __exportAll } from "../_runtime.mjs";
import {
  n as require_jsx_runtime,
  t as QueryClientProvider,
} from "../_libs/react+tanstack__react-query.mjs";
import {
  c as HeadContent,
  d as Outlet,
  f as lazyRouteComponent,
  g as useRouter,
  h as Link,
  m as createRootRouteWithContext,
  p as createFileRoute,
  s as Scripts,
  u as createRouter,
} from "../_libs/@tanstack/react-router+[...].mjs";
import { t as QueryClient } from "../_libs/tanstack__query-core.mjs";
//#region node_modules/.nitro/vite/services/ssr/assets/router-DkXKgi1X.js
var router_DkXKgi1X_exports = /* @__PURE__ */ __exportAll({
  getRouter: () => getRouter,
});
var import_jsx_runtime = require_jsx_runtime();
var styles_default = "/assets/styles-B25xTD3Z.css";
function NotFoundComponent() {
  return /* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
    className:
      "flex min-h-screen items-center justify-center bg-background px-4",
    children: /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
      className: "max-w-md text-center",
      children: [
        /* @__PURE__ */ (0, import_jsx_runtime.jsx)("h1", {
          className: "text-7xl font-bold text-foreground",
          children: "404",
        }),
        /* @__PURE__ */ (0, import_jsx_runtime.jsx)("h2", {
          className: "mt-4 text-xl font-semibold text-foreground",
          children: "Page not found",
        }),
        /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
          className: "mt-2 text-sm text-muted-foreground",
          children:
            "The page you're looking for doesn't exist or has been moved.",
        }),
        /* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
          className: "mt-6",
          children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Link, {
            to: "/",
            className:
              "inline-flex items-center justify-center rounded-md bg-primary px-4 py-2 text-sm font-medium text-primary-foreground transition-colors hover:bg-primary/90",
            children: "Go home",
          }),
        }),
      ],
    }),
  });
}
function ErrorComponent({ error, reset }) {
  console.error(error);
  const router = useRouter();
  return /* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
    className:
      "flex min-h-screen items-center justify-center bg-background px-4",
    children: /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
      className: "max-w-md text-center",
      children: [
        /* @__PURE__ */ (0, import_jsx_runtime.jsx)("h1", {
          className: "text-xl font-semibold tracking-tight text-foreground",
          children: "This page didn't load",
        }),
        /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
          className: "mt-2 text-sm text-muted-foreground",
          children:
            "Something went wrong on our end. You can try refreshing or head back home.",
        }),
        /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
          className: "mt-6 flex flex-wrap justify-center gap-2",
          children: [
            /* @__PURE__ */ (0, import_jsx_runtime.jsx)("button", {
              onClick: () => {
                router.invalidate();
                reset();
              },
              className:
                "inline-flex items-center justify-center rounded-md bg-primary px-4 py-2 text-sm font-medium text-primary-foreground transition-colors hover:bg-primary/90",
              children: "Try again",
            }),
            /* @__PURE__ */ (0, import_jsx_runtime.jsx)("a", {
              href: "/",
              className:
                "inline-flex items-center justify-center rounded-md border border-input bg-background px-4 py-2 text-sm font-medium text-foreground transition-colors hover:bg-accent",
              children: "Go home",
            }),
          ],
        }),
      ],
    }),
  });
}
var Route$1 = createRootRouteWithContext()({
  head: () => ({
    meta: [
      { charSet: "utf-8" },
      {
        name: "viewport",
        content: "width=device-width, initial-scale=1",
      },
      { title: "DissolveAI — Molecular Solubility Prediction" },
      {
        name: "description",
        content:
          "Machine-learning prediction of molecular solubility across solvents using RDKit features and a Random Forest regressor.",
      },
      {
        name: "author",
        content: "Iram Javed",
      },
      {
        property: "og:title",
        content: "DissolveAI — Molecular Solubility Prediction",
      },
      {
        property: "og:description",
        content:
          "Machine-learning prediction of molecular solubility across solvents using RDKit features and a Random Forest regressor.",
      },
      {
        property: "og:type",
        content: "website",
      },
      {
        name: "twitter:card",
        content: "summary_large_image",
      },
    ],
    links: [
      {
        rel: "stylesheet",
        href: styles_default,
      },
      {
        rel: "preconnect",
        href: "https://fonts.googleapis.com",
      },
      {
        rel: "preconnect",
        href: "https://fonts.gstatic.com",
        crossOrigin: "anonymous",
      },
      {
        rel: "stylesheet",
        href: "https://fonts.googleapis.com/css2?family=Newsreader:ital,opsz,wght@0,6..72,300;0,400;0,500;1,6..72,400&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap",
      },
      {
        rel: "icon",
        href: "/favicon.ico",
        type: "image/x-icon",
      },
    ],
  }),
  shellComponent: RootShell,
  component: RootComponent,
  notFoundComponent: NotFoundComponent,
  errorComponent: ErrorComponent,
});
function RootShell({ children }) {
  return /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("html", {
    lang: "en",
    children: [
      /* @__PURE__ */ (0, import_jsx_runtime.jsx)("head", {
        children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)(HeadContent, {}),
      }),
      /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("body", {
        children: [
          children,
          /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Scripts, {}),
        ],
      }),
    ],
  });
}
function RootComponent() {
  const { queryClient } = Route$1.useRouteContext();
  return /* @__PURE__ */ (0, import_jsx_runtime.jsx)(QueryClientProvider, {
    client: queryClient,
    children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Outlet, {}),
  });
}
var $$splitComponentImporter = () => import("./routes-BOwAXmoD.mjs");
var rootRouteChildren = {
  IndexRoute: createFileRoute("/")({
    head: () => ({
      meta: [
        { title: "DissolveAI — Molecular Solubility Prediction" },
        {
          name: "description",
          content:
            "DissolveAI estimates compound solubility across 70 solvents from SMILES structures using RDKit features and a Random Forest regressor.",
        },
        {
          property: "og:title",
          content: "DissolveAI — Molecular Solubility Prediction",
        },
        {
          property: "og:description",
          content:
            "Predict molecular solubility before the experiment: RDKit descriptors, Morgan fingerprints, solvent encoding, Random Forest regression.",
        },
      ],
    }),
    component: lazyRouteComponent($$splitComponentImporter, "component"),
  }).update({
    id: "/",
    path: "/",
    getParentRoute: () => Route$1,
  }),
};
var routeTree = Route$1._addFileChildren(rootRouteChildren)._addFileTypes();
var getRouter = () => {
  const queryClient = new QueryClient();
  return createRouter({
    routeTree,
    context: { queryClient },
    scrollRestoration: true,
    defaultPreloadStaleTime: 0,
  });
};
//#endregion
export { getRouter, router_DkXKgi1X_exports as t };
